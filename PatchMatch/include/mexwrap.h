/* 
 * File:   mexwrap.h
 * Author: xion
 *
 * Created on April 28, 2014, 12:31 AM
 */

#ifndef MEXWRAP_H
#define	MEXWRAP_H

#include "mexutil.h"
#include "algebra.h"

namespace pm {

	/**
	* Image matrix wrapper around MxN matlab matrices
	*/
	struct MatWrapper {
	   union {
		   int rows;
		   int height;
	   };
	   union {
		   int cols;
		   int width;
	   };
	   //! returns the type of the matrix elements
	   inline int type() const {
		   return flags;
	   }
	   //! returns the depth of the matrix elements
	   inline int depth() const {
		   return IM_MAT_DEPTH(flags);
	   }
	   //! returns the number of channels in each matrix element
	   inline int channels() const {
		   return IM_MAT_CN(flags);
	   }

	   MatWrapper() : flags(IM_UNKNOWN), data(NULL), mx(NULL){
	   }

	   MatWrapper(int h, int w, int dataType) : height(h), width(w), flags(dataType){
		   int elemSize = IM_SIZEOF_BY_CHANNEL(dataType);
		   int num_ch = IM_MAT_CN(dataType);
		   mx = num_ch > 1 ? mxCreateMatrix(h, w, num_ch, classOf(dataType))
				   : mxCreateMatrix(h, w, classOf(dataType));
		   data = reinterpret_cast<byte *>(mxGetData(mx));
		   if(data != NULL){
			   step[0] = elemSize;
			   step[1] = height * elemSize;
			   step[2] = height * width * elemSize;
		   } else {
			   std::cerr << "No byte for dataType=" << dataType << " with size=" << elemSize << "\n";
			   step[0] = step[1] = step[2] = 0;
		   }
	   }
	   MatWrapper(int h, int w, mxClassID c, int num_ch = 1) : height(h), width(w) {
		   mx = num_ch > 1 ? mxCreateMatrix(h, w, num_ch, c)
				   : mxCreateMatrix(h, w, c);
		   data = reinterpret_cast<byte *>(mxGetData(mx));
		   flags = IM_MAKETYPE(depthOf(c), num_ch);
		   int elemSize = IM_SIZEOF_BY_CHANNEL(flags);
		   if(data != NULL){
			   step[0] = elemSize;
			   step[1] = height * elemSize;
			   step[2] = height * width * elemSize;
		   }
	   }
	   
	   MatWrapper(const mxArray *arr)
	   : height(mxGetDimensions(arr)[0]), width(mxGetDimensions(arr)[1]) {
		   int num_ch = mxGetNumberOfDimensions(arr) < 3 ? 1 : mxGetDimensions(arr)[2];
		   flags = IM_MAKETYPE(depthOf(arr), num_ch);
		   mx = const_cast<mxArray *>(arr); // you shouldn't try modifying this instance!
		   data = reinterpret_cast<byte *>(mxGetData(mx));
		   int elemSize = IM_SIZEOF_BY_CHANNEL(flags);
		   if(data != NULL){
			   step[0] = elemSize;
			   step[1] = height * elemSize;
			   step[2] = height * width * elemSize;
		   } else {
			   std::cerr << "No byte for dataType=" << flags << " with size=" << elemSize << "\n";
			   step[0] = step[1] = step[2] = 0;
		   }
	   }
	   operator mxArray*() {
		   return mx;
	   }
	   operator const mxArray*() const {
		   return mx;
	   }

	   inline bool empty() const {
		   return !data;
	   }

	   //! Pointer access
	   template <typename T>
	   inline const T *ptr(int y, int x, int ch = 0) const {
		   // std::cout << "ptr(" << y << ", " << x << ") -> x [" << step[1] << ", " << step[0] << "] const\n";  
		   return reinterpret_cast<const T*>(data + y * step[0] + x * step[1] + ch * step[2]);
	   }
	   template <typename T>
	   inline T *ptr(int y, int x, int ch = 0) {
		   // std::cout << "ptr(" << y << ", " << x << ") -> x [" << step[1] << ", " << step[0] << "]\n";
		   return reinterpret_cast<T*>(data + y * step[0] + x * step[1] + ch * step[2]);
	   }

	   //! Element access
	   template <typename T>
	   inline const T &at(int y, int x, int ch = 0) const {
		   if(x < 0 || x >= width || y < 0 || y >= height || ch < 0 || ch >= channels()) {
				std::cerr << "@" << y << "/" << x << "/" << ch << " of ";
				std::cerr << height << "/" << width << "/" << channels() << "\n";
				mexErrMsgIdAndTxt("MATLAB:img:ptr", "Ref out of image bounds!");
		   }
		   return *ptr<T>(y, x, ch);
	   }
	   template <typename T>
	   inline T &at(int y, int x, int ch = 0) {
		   if(x < 0 || x >= width || y < 0 || y >= height || ch < 0 || ch >= channels()) {
				std::cerr << "@" << y << "/" << x << "/" << ch << " of ";
				std::cerr << height << "/" << width << "/" << channels() << "\n";
				mexErrMsgIdAndTxt("MATLAB:img:ptr", "Ref out of image bounds!");
		   }
		   return *ptr<T>(y, x, ch);
	   }
	   
	   typedef unsigned char uchar;
	   typedef signed char schar;
	   typedef unsigned int uint;
	   typedef signed int sint;
	   typedef unsigned long ulong;
	   typedef signed long slong;
	   
	   template <typename T>
	   inline const T read(int y, int x, int ch = 0) const {
		   if(x < 0 || x >= width || y < 0 || y >= height || ch < 0 || ch >= channels()) {
				std::cerr << "r@" << y << "/" << x << "/" << ch << " of ";
				std::cerr << height << "/" << width << "/" << channels() << "\n";
				mexErrMsgIdAndTxt("MATLAB:img:ptr", "Ref out of image bounds!");
		   }
		   switch(mxGetClassID(mx)) {
			   case mxSINGLE_CLASS: return T(*ptr<float>(y, x, ch));
			   case mxDOUBLE_CLASS: return T(*ptr<double>(y, x, ch));
			   case mxUINT8_CLASS: return T(*ptr<uchar>(y, x, ch));
			   case mxINT8_CLASS: return T(*ptr<schar>(y, x, ch));
			   case mxUINT32_CLASS: return T(*ptr<uint>(y, x, ch));
			   case mxINT32_CLASS: return T(*ptr<sint>(y, x, ch));
			   case mxUINT64_CLASS: return T(*ptr<ulong>(y, x, ch));
			   case mxINT64_CLASS: return T(*ptr<slong>(y, x, ch));
			   case mxLOGICAL_CLASS: return T(*ptr<bool>(y, x, ch));
			   default:
				   mexErrMsgIdAndTxt("MATLAT:img:ptr", "Unsupported class %d", mxGetClassID(mx));
				   return T(0);
		   }
	   }
	   template <typename T>
	   inline void update(int y, int x, int ch, T value) {
		   if(x < 0 || x >= width || y < 0 || y >= height || ch < 0 || ch >= channels()) {
				std::cerr << "u@" << y << "/" << x << "/" << ch << " of ";
				std::cerr << height << "/" << width << "/" << channels() << "\n";
				mexErrMsgIdAndTxt("MATLAB:img:ptr", "Ref out of image bounds!");
		   }
		   switch(mxGetClassID(mx)) {
			   case mxSINGLE_CLASS: *ptr<float>(y, x, ch) = float(value); break;
			   case mxDOUBLE_CLASS: *ptr<double>(y, x, ch) = double(value); break;
			   case mxUINT8_CLASS: *ptr<uchar>(y, x, ch) = uchar(value); break;
			   case mxINT8_CLASS: *ptr<schar>(y, x, ch) = schar(value); break;
			   case mxUINT32_CLASS: *ptr<uint>(y, x, ch) = uint(value); break;
			   case mxINT32_CLASS: *ptr<sint>(y, x, ch) = sint(value); break;
			   case mxUINT64_CLASS: *ptr<ulong>(y, x, ch) = ulong(value); break;
			   case mxINT64_CLASS: *ptr<slong>(y, x, ch) = slong(value); break;
			   case mxLOGICAL_CLASS: *ptr<bool>(y, x, ch) = bool(value); break;
			   default:
				   mexErrMsgIdAndTxt("MATLAT:img:ptr", "Unsupported class %d", mxGetClassID(mx));
		   }
	   }
	   template <typename T>
	   inline void update(int y, int x, T value) {
		   update(y, x, 0, value);
	   }

	private:
	   int flags;
	   byte *data;
	   mxArray *mx;
	   int step[3];
	};
	
	typedef MatWrapper MatXD;

}

#endif	/* MEXWRAP_H */

