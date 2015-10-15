/*******************************************************************************
 * mexutil.h - mex helpers
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef MEXUTIL_H
#define	MEXUTIL_H

// omp
#if _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_set_num_threads(x)
#endif
// other
#include <cstring>
#include <mex.h>
#include <filter.h>
#include <gpm.h>
#include <nnf.h>
#include <vector>
#include <voting/histogram.h>

using pm::Image;
using pm::DataDepth;

template <typename Scalar>
mxClassID classID();

inline mxClassID classOf(int dataType) {
	int depth = IM_MAT_DEPTH(dataType);
	switch (depth) {
		case IM_8U: return mxUINT8_CLASS;
		case IM_8S: return mxINT8_CLASS;
		case IM_32S: return mxINT32_CLASS;
		case IM_32F: return mxSINGLE_CLASS;
		case IM_64F: return mxDOUBLE_CLASS;
		default:
			mexErrMsgIdAndTxt("MATLAB:mex:classOf", "Unsupported class!");
			return mxUNKNOWN_CLASS;
	}
}
inline mxClassID classOf(const Image &img) {
	return classOf(img.depth());
}

inline int depthOf(mxClassID c) {
	switch(c) {
		case mxUINT8_CLASS: return IM_8U;
		case mxINT8_CLASS: return IM_8S;
		case mxINT32_CLASS: return IM_32S;
		case mxSINGLE_CLASS: return IM_32F;
		case mxDOUBLE_CLASS: return IM_64F;
		case mxLOGICAL_CLASS: return IM_8U;
		default:
			mexErrMsgIdAndTxt("MATLAB:mex:depthOf", "Unsupported depth!");
			return -1;
	}
}
inline int depthOf(const mxArray *arr) {
	return depthOf(mxGetClassID(arr));
}

inline bool mxStringEquals(const mxArray *A, const char *s) {
	char buf[256];
	if (!mxIsChar(A)) {
		return false;
	}
	if (mxGetString(A, buf, 255)) {
		return false;
	}
	return strcmp(s, buf) == 0;
}

inline mxArray *mxCreateMatrix(int rows, int cols, mxClassID type = mxSINGLE_CLASS) {
	mwSize sz[2] = {rows, cols};
	return mxCreateNumericArray(2, sz, type, mxREAL);
}

inline mxArray *mxCreateMatrix(int rows, int cols, int channels, mxClassID type = mxSINGLE_CLASS) {
	mwSize sz[3] = {rows, cols, channels};
	return mxCreateNumericArray(3, sz, type, mxREAL);
}

template <typename Scalar>
inline mxArray *mxCreateMatrix(const Image &img) {
	return mxCreateMatrix(img.rows, img.cols, img.channels(), classID<Scalar>());
}

inline mxArray *mxCreateMatrix(const Image &img) {
	switch (img.depth()) {
		case IM_8U: return mxCreateMatrix(img.rows, img.cols, img.channels(), mxUINT8_CLASS);
		case IM_32F: return mxCreateMatrix(img.rows, img.cols, img.channels(), mxSINGLE_CLASS);
		case IM_64F: return mxCreateMatrix(img.rows, img.cols, img.channels(), mxDOUBLE_CLASS);
		default:
			mexErrMsgIdAndTxt("MATLAB:mex:createMatrixFor", "Class type not supported!");
			return NULL;
	}
}

template <typename Scalar>
inline mxArray *mxCreateScalar(Scalar s) {
	mxArray *d = mxCreateNumericArray(1, 1, classID<Scalar>());
	Scalar *ptr = reinterpret_cast<Scalar *>(mxGetData(d));
	*ptr = s;
	return d;
}

inline double mxCheckedScalar(const mxArray *a, const char *s) {
	if (!mxIsNumeric(a)) {
		std::cerr << "ClassID = " << mxGetClassID(a) << "\n";
		mexErrMsgIdAndTxt("MATLAB:mex:invalidInput", s);
	}
	return mxGetScalar(a);
}

inline bool mxHasField(const mxArray *options, mwIndex i, const char *fieldname) {
	if(!mxIsStruct(options)){
        mexErrMsgIdAndTxt("MATLAB:mex:mxHasField", "Invalid call mxHasField for '%s' on non-struct object.", fieldname);
        return false;
    }
	return mxGetField(options, i, fieldname) != NULL;
}

inline double mxScalarField(const mxArray *options, mwIndex i, const char *fieldname, double defaultValue = 0) {
	if(!mxIsStruct(options)){
        mexErrMsgIdAndTxt("MATLAB:mex:mxHasField", "Invalid call mxHasField for '%s' on non-struct object.", fieldname);
        return defaultValue;
    }
	const mxArray *tmp = mxGetField(options, i, fieldname);
	if(tmp == NULL)
		return defaultValue;
	if(!mxIsNumeric(tmp))
		return defaultValue;
	return mxGetScalar(tmp);
}

inline bool mxBoolField(const mxArray *options, mwIndex i, const char *fieldname, bool defValue = false) {
    if(!mxIsStruct(options)){
        mexErrMsgIdAndTxt("MATLAB:mex:mxBoolField", "Invalid call mxBoolField for '%s' on non-struct object.", fieldname);
        return defValue;
    }
	const mxArray *tmp = mxGetField(options, i, fieldname);
	if(!tmp) return defValue;
	if(mxIsNumeric(tmp)){
		return mxGetScalar(tmp) != 0;
	} else if(mxIsLogical(tmp) && mxGetNumberOfElements(tmp) > 0) {
		return mxGetLogicals(tmp)[0];
	} else {
		mexErrMsgIdAndTxt("MATLAB:mex:invalidInput", "Parameter %s is not a valid boolean expression.", fieldname);
		return defValue;
	}
}

template <typename Scalar>
inline void mxLoadScalars(Scalar *ptr, int n, const mxArray *arr, const char *s) {
	if (!mxIsNumeric(arr)) {
		mexErrMsgIdAndTxt("MATLAB:mex:invalidInput", s);
	}
	if (mxGetClassID(arr) != classID<Scalar>()) {
		mexErrMsgIdAndTxt("MATLAB:mex:invalidInput", s);
	}
	const Scalar *data = reinterpret_cast<const Scalar *> (mxGetData(arr));
	switch (mxGetNumberOfElements(arr)) {
		case 1:
			std::fill(ptr, ptr + n, Scalar(mxGetScalar(arr)));
			break;
		case 5:
			for (int i = 0; i < 5; ++i) ptr[i] = data[i];
			break;
		default:
			mexErrMsgIdAndTxt("MATLAB:mex:invalidInput", s);
			break;
	}
}

inline void mxLoadFilter(pm::Filter &filter, const mxArray *arr, const char *s) {
		if (!mxIsNumeric(arr)) {
			mexErrMsgIdAndTxt("MATLAB:vote:wrongWeights",
					"vote_weight should be an array of numbers!");
		}
		int m = mxGetM(arr);
		int n = mxGetN(arr);
        if (m == n && m == 1) {
            filter = pm::Filter(filter.width, mxGetScalar(arr));
            return;
        }
		if (m != n || m != filter.width) {
			mexErrMsgIdAndTxt("MATLAB:filter:wrongSize", s);
		}
		switch (mxGetClassID(arr)) {
			case mxSINGLE_CLASS:
			{
				float *data = reinterpret_cast<float *> (mxGetData(arr));
				for (int y = 0; y < m; ++y) {
					for (int x = 0; x < n; ++x) {
						filter[y][x] = data[y + x * m];
					}
				}
			}
				break;
			case mxDOUBLE_CLASS:
			{
				double *data = reinterpret_cast<double *> (mxGetData(arr));
				for (int y = 0; y < m; ++y) {
					for (int x = 0; x < n; ++x) {
						filter[y][x] = float(data[y + x * m]);
					}
				}
			}
				break;
			default:
				mexErrMsgIdAndTxt("MATLAB:vote:voteClass",
						"vote_weight must be single or double.");
				break;
		}
}

template <typename S, typename T>
inline void mxLoadVector(std::vector<T> &v, const mxArray *arr, const char *s) {
	if (!mxIsNumeric(arr)) {
		mexErrMsgIdAndTxt("MATLAB:mex:invalidInput", s);
	}
	const S *data = reinterpret_cast<const S *> (mxGetData(arr));
	int N = mxGetNumberOfElements(arr);
	if(N == 0) mexErrMsgIdAndTxt("MATLAB:mex:mxLoadVector", s);
	// resize the vector and fill it
	v.resize(N);
	for(int i = 0; i < N; ++i) v[i] = T(data[i]);
}

template <typename T>
inline void mxLoadVector(std::vector<T> &v, const mxArray *arr, const char *s) {
	switch (mxGetClassID(arr)) {
		case mxINT8_CLASS: mxLoadVector<char, T>(v, arr, s); break;
		case mxUINT8_CLASS: mxLoadVector<unsigned char, T>(v, arr, s); break;
		case mxSINGLE_CLASS: mxLoadVector<float, T>(v, arr, s); break;
		case mxDOUBLE_CLASS: mxLoadVector<double, T>(v, arr, s); break;
		default:
			mexErrMsgIdAndTxt("MATLAB:mex:mxLoadVector", "Class type not supported!");
	}
}

template <typename GB>
inline void mxLoadGB(const mxArray *options) {
	typedef typename GB::Vec3 Vec3;
	mxArray *tmp;
	if (tmp = mxGetField(options, 0, "min_bias")) {
		mxLoadScalars(GB::minBias.data, 3, tmp, "Invalid min_bias!");
	} else {
		GB::minBias = Vec3(-10, -50, -50);
	}
	if (tmp = mxGetField(options, 0, "max_bias")) {
		mxLoadScalars(GB::maxBias.data, 3, tmp, "Invalid max_bias!");
	} else {
		GB::maxBias = Vec3(10, 50, 50);
	}
	if (tmp = mxGetField(options, 0, "min_bias")) {
		mxLoadScalars(GB::minGain.data, 3, tmp, "Invalid min_gain!");
	} else {
		GB::minGain = Vec3(0.9, 0.9, 0.9);
	}
	if (tmp = mxGetField(options, 0, "max_gain")) {
		mxLoadScalars(GB::maxGain.data, 3, tmp, "Invalid max_gain!");
	} else {
		GB::maxGain = Vec3(1.5, 1.5, 1.5);
	}
}

template <typename Scalar>
inline mxArray *mxImageToArray(const Image &img) {
	// std::cout << "mxImageToArray\n";
	mxArray *arr = mxCreateMatrix<Scalar>(img);
	const int offset = img.rows * img.cols;
	Scalar *data = reinterpret_cast<Scalar *> (mxGetData(arr));
	// #pragma omp parallel for collapse(2)
	for (int y = 0; y < img.rows; ++y) {
		for (int x = 0; x < img.cols; ++x) {
			const Scalar *iptr = img.ptr<Scalar>(y, x);
			for (int ch = 0; ch < img.channels(); ++ch) {
				// transposing!
				data[y + x * img.rows + offset * ch] = iptr[ch];
			}
		}
	}
	return arr;
	// return MxArray(img, mxUNKNOWN_CLASS, false);
}

inline mxArray *mxImageToArray(const Image &img) {
	switch (img.depth()) {
		case IM_8S: return mxImageToArray<char>(img);
		case IM_8U: return mxImageToArray<unsigned char>(img);
		case IM_32F: return mxImageToArray<float>(img);
		case IM_64F: return mxImageToArray<double>(img);
		default:
			mexErrMsgIdAndTxt("MATLAB:mex:mxImageToArray", "Class type not supported!");
			return NULL;
	}
}

template <typename Scalar>
inline void mxCheckImage(const Image &img, const char *errMsg = "Corrupted image with pixels out of bounds!") {
	bool err = false;
	for (int y = 0; y < img.rows; ++y) {
		for (int x = 0; x < img.cols; ++x) {
			const Scalar *ptr = img.ptr<Scalar>(y, x);
			for (int c = 0; c < img.channels(); ++c) {
				Scalar d = ptr[c];
				if ((d + 1) == d) {
					std::cout << "Invalid p@" << y << "/" << x << "/c" << c << " = " << d << "\n";
					err = true;
				}
			}
		}
	}
	if (err) mexErrMsgIdAndTxt("MATLAB:mex:invalid_image", errMsg);
}

template <typename Scalar>
inline Image mxArrayToImage(const mxArray *arr, const char *errMsg) {
	if (!mxIsNumeric(arr)) {
		mexErrMsgIdAndTxt("MATLAB:mex:invalidInput", "Invalid image array.");
	}
	int h = mxGetDimensions(arr)[0];
	int w = mxGetDimensions(arr)[1];
	int num_ch = mxGetNumberOfDimensions(arr) < 3 ? 1 : mxGetDimensions(arr)[2];
	const int offset = h * w; 
	// std::cout << "h=" << h << ", w=" << w << ", ch=" << num_ch << ", offset=" << offset << "\n";
	assert(mxGetClassID(arr) == classID<Scalar>());
	Image img(h, w, IM_MAKETYPE(DataDepth<Scalar>::value, num_ch));
	// std::cout << "depth=" << DataDepth<Scalar>::value << ", ch=" << num_ch << "=" << img.channels() << "\n";
	const Scalar *data = reinterpret_cast<const Scalar *> (mxGetData(arr));
	// std::cout << "steps: " << img.step[0] << ", " << img.step[1] << ", " << img.step[2] << "\n";
	// /!\ img.step[2] might be wrong! do not use at(y, x, ch) indexing!
	// => use img.ptr(y, x)
	if (num_ch == 1) {
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < img.rows; ++y) {
			for (int x = 0; x < img.cols; ++x) {
				// transposing!
				img.at<Scalar>(y, x) = data[y + x * img.rows];
			}
		}
	} else {
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < img.rows; ++y) {
			for (int x = 0; x < img.cols; ++x) {
				Scalar *iptr = img.ptr<Scalar>(y, x);
				for (int ch = 0; ch < img.channels(); ++ch) {
					// std::cout << "at(" << y << ", " << x << ", " << ch << ") = ";
					// transposing!
					Scalar v = data[y + x * img.rows + offset * ch];
					// std::cout << v << "\n";
					iptr[ch] = v;
				}
			}
		}
	}
	mxCheckImage<Scalar>(img, errMsg);
	return img;
	// return MxArray(arr).toMat(CV_USRTYPE1, false);
}

inline Image mxArrayToImage(const mxArray *arr, const char *err = "Corrupted image with pixels out of bounds!") {
	switch (mxGetClassID(arr)) {
		case mxINT8_CLASS: return mxArrayToImage<char>(arr, err);
		case mxUINT8_CLASS: return mxArrayToImage<unsigned char>(arr, err);
		case mxSINGLE_CLASS: return mxArrayToImage<float>(arr, err);
		case mxDOUBLE_CLASS: return mxArrayToImage<double>(arr, err);
		default:
			mexErrMsgIdAndTxt("MATLAB:mex:mxArrayToImage", "Class type not supported!");
			return Image();
	}
}

inline void mxLoadTexture(pm::Texture &tex, const mxArray *arr, const char *err = "Corrupted texture!") {
	if(mxIsCell(arr)) {
		int N = mxGetNumberOfElements(arr);
		if (N < 1) {
			mexErrMsgIdAndTxt("MATLAB:mex:mxLoadTexture", "Empty texture cell!");
		}
		const mxArray *cell = mxGetCell(arr, 0);
		if(mxIsNumeric(cell))
			tex = pm::Texture(mxArrayToImage(cell, "Invalid first texture cell!"), N - 1); //< stack is allocated here!
		else
			mexErrMsgIdAndTxt("MATLAB:mex:mxLoadTexture", "Texture cell is not numeric!");
		for(int i = 1; i < N; ++i) {
			cell = mxGetCell(arr, i);
			if(cell == NULL){
				mexErrMsgIdAndTxt("MATLAB:mex:mxLoadTexture", "Null cell!");
			}
			tex.stack[i - 1] = mxArrayToImage(cell);
			// check image state
			// std::cout << "T" << i << ": "  << tex.stack[i-1].at<float>(tex.height - 1, tex.width - 1) << "\n";
		}
		// std::cout << "loaded texture (d=" << (N-1) << ")" << "\n";
	} else if(mxIsNumeric(arr)) {
		tex = pm::Texture(mxArrayToImage(arr, err));
	} else {
		mexErrMsgIdAndTxt("MATLAB:mex:mxLoadTexture", "Unknown texture class!");
	}
}

template <typename Patch, typename Scalar>
inline void mxLoadNNF(pm::NearestNeighborField<Patch, Scalar> &nnf, const mxArray *narr,
		const char *numericErrorMsg = "The NNF should contain numbers!") {
	if (!mxIsNumeric(narr)) {
		mexErrMsgIdAndTxt("MATLAB:vote:invalidArgument",
				numericErrorMsg);
	}
	const int h = mxGetDimensions(narr)[0];
	const int w = mxGetDimensions(narr)[1];
	const int num_ch = mxGetNumberOfDimensions(narr) < 3 ? 1 : mxGetDimensions(narr)[2];
	const int offset = h * w;
	// check that sizes match
	if (h != nnf.height || w != nnf.width) {
		mexErrMsgIdAndTxt("MATLAB:mex:nnf_size",
				"Invalid nnf size, got HxW=%dx%d instead of %dx%d", h, w, nnf.height, nnf.width);
	}
	// check that the class is valid
	if (mxGetClassID(narr) != classID<float>()) {
		mexErrMsgIdAndTxt("MATLAB:mex:nnf_class",
				"Invalid NNF class, should be single!");
	}
	// check the channels
	if (num_ch < 2) {
		mexErrMsgIdAndTxt("MATLAB:mex:nnf_channels",
				"The NNF should have at least two channels!");
	}
	// load the real nnf data
	const float *fptr = reinterpret_cast<const float *> (mxGetData(narr));
#if _OPENMP
#pragma omp parallel for collapse(2) ordered
#endif
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			const float *base = fptr + y + x * h; // transposed!
//			std::cout << "f@" << y << "/" << x << ": ";
//			for(int i = 0; i < num_ch; ++i) std::cout << base[i * offset] << ", ";
//			std::cout << "\n";
			Patch &patch = nnf.get(y, x);
			patch.load(base, num_ch, offset);
//			std::cout << "- num_ch=" << num_ch << "\n";
//			std::cout << "p: " << patch.y << "/" << patch.x << "\n";
//			float d[10];
//			float *data = &d[0];
//			patch.store(data, num_ch);
//			std::cout << "p@" << y << "/" << x << ": ";
//			for(int i = 0; i < num_ch; ++i) std::cout << d[i] << ", ";
//			std::cout << "\n";
			
			//				if((x + y) % 100 == 0){
			//				cv::Vec<float, 5> d;
			//				patch.store(d, 5);
			//				std::cout << "@" << y << "/" << x << ": " << d << "\n";
			//				}
			//^ patch done
		}
	}
}

// boring implementations

template <>
mxClassID classID<float>() {
	return mxSINGLE_CLASS;
}

template <>
mxClassID classID<double>() {
	return mxDOUBLE_CLASS;
}

template <>
mxClassID classID<unsigned char>() {
	return mxUINT8_CLASS;
}

template <>
mxClassID classID<char>() {
	return mxINT8_CLASS;
}

#endif	/* MEXUTIL_H */

