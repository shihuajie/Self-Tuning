/*******************************************************************************
 * algebra.h - matrix and vectors for imaging
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef ALGEBRA_H
#define	ALGEBRA_H

#include <boost/shared_ptr.hpp>
#include <cmath>
#include <iostream>
#include <mex.h>

// bit boundaries
#define IM_CN_MAX     512
#define IM_CN_SHIFT   3
#define IM_DEPTH_MAX  (1 << IM_CN_SHIFT)

// bits for the type
#define IM_8U   0
#define IM_8S   1
#define IM_32S  2
#define IM_32F  3
#define IM_64F  4
#define IM_USRTYPE 7
#define IM_UNKNOWN -1
// ... we can add two more types

// depth of type
#define IM_MAT_DEPTH_MASK       (IM_DEPTH_MAX - 1)
#define IM_MAT_DEPTH(flags)     ((flags) & IM_MAT_DEPTH_MASK)

// full type generation
#define IM_MAKETYPE(depth, cn) (IM_MAT_DEPTH(depth) + (((cn)-1) << IM_CN_SHIFT))
#define IM_MAKE_TYPE(depth, cn) IM_MAKETYPE(depth, cn)

#define IM_8UC1 IM_MAKETYPE(IM_8U,1)
#define IM_8UC3 IM_MAKETYPE(IM_8U,3)
#define IM_8UC(n) IM_MAKETYPE(IM_8U,(n))

#define IM_32SC1 IM_MAKETYPE(IM_32S,1)
#define IM_32SC3 IM_MAKETYPE(IM_32S,3)
#define IM_32SC(n) IM_MAKETYPE(IM_32S,(n))

#define IM_32FC1 IM_MAKETYPE(IM_32F,1)
#define IM_32FC3 IM_MAKETYPE(IM_32F,3)
#define IM_32FC(n) IM_MAKETYPE(IM_32F,(n))

#define IM_64FC1 IM_MAKETYPE(IM_64F,1)
#define IM_64FC3 IM_MAKETYPE(IM_64F,3)
#define IM_64FC(n) IM_MAKETYPE(IM_64F,(n))

// channel & type extraction
#define IM_MAT_CN_MASK          ((IM_CN_MAX - 1) << IM_CN_SHIFT)
#define IM_MAT_CN(flags)        ((((flags) & IM_MAT_CN_MASK) >> IM_CN_SHIFT) + 1)
#define IM_MAT_TYPE_MASK        (IM_DEPTH_MAX*IM_CN_MAX - 1)
#define IM_MAT_TYPE(flags)      ((flags) & IM_MAT_TYPE_MASK)

// size of a type
#define IM_SIZEOF_DEPTH(depth)	(depth <= IM_8S ? sizeof(byte) : \
									(depth <= IM_32F ? sizeof(int) : \
									(depth == IM_64F ? sizeof(double) : sizeof(int*)) ))
#define IM_SIZEOF_IMPL(depth, cn)	(cn > 0 ? cn * IM_SIZEOF_DEPTH(depth) : 0)
#define IM_SIZEOF(flags)		IM_SIZEOF_IMPL(IM_MAT_DEPTH(flags), IM_MAT_CN(flags))
#define IM_SIZEOF_BY_CHANNEL(flags)	IM_SIZEOF_IMPL(IM_MAT_DEPTH(flags), 1)

namespace pm {
	
	/**
	 * Basic byte type
	 */
	typedef unsigned char byte;

	/**
	 * Single data type numerical representation
	 */
	template<typename Scalar> struct DataDepth {};

	template<> struct DataDepth<bool> { enum { value = IM_8U, fmt=(int)'u', exact = 1 }; };
	template<> struct DataDepth<unsigned char> { enum { value = IM_8U, fmt=(int)'u', exact = 1 }; };
	template<> struct DataDepth<signed char> { enum { value = IM_8S, fmt=(int)'c', exact = 1 }; };
	template<> struct DataDepth<char> { enum { value = IM_8S, fmt=(int)'c', exact = 1 }; };
	template<> struct DataDepth<int> { enum { value = IM_32S, fmt=(int)'i', exact = 1 }; };
	template<> struct DataDepth<unsigned int> { enum { value = IM_32S, fmt=(int)'i', exact = 1 }; };
	template<> struct DataDepth<float> { enum { value = IM_32F, fmt=(int)'f', exact = 0 }; };
	template<> struct DataDepth<double> { enum { value = IM_64F, fmt=(int)'d', exact = 0 }; };
	template<typename Scalar> struct DataDepth<Scalar*> { enum { value = IM_USRTYPE, fmt=(int)'r', exact = 1 }; };
	
	template<int flag> struct DataSize {};

	/**
	 * \brief Small vector that stays on the stack
	 */
	template<typename T, int cn>
	struct Vec {
		typedef T scalar;
		typedef Vec<T, cn> vec;

		enum {
			typeDepth = DataDepth<scalar>::value,
			channels = cn,
			type = IM_MAKETYPE(typeDepth, channels)
		};

		//! default constructor

		Vec() {
			// note: we do not initialize data here!
		}

		Vec(T v0) {
			data[0] = v0;
		}

		Vec(T v0, T v1) {
			data[0] = v0;
			data[1] = v1;
		}

		Vec(T v0, T v1, T v2) {
			data[0] = v0;
			data[1] = v1;
			data[2] = v2;
		}

		explicit Vec(const T* vals) {
			for (int i = 0; i < channels; ++i) data[i] = vals[i];
		}

		static vec all(T alpha) {
			vec m;
			std::fill(m.data, m.data + channels, alpha);
			return m;
		}

		inline static vec zeros() {
			return all(0);
		}

		inline static vec ones() {
			return all(1);
		}

		//! convertion to another data type

		template<typename T2> operator Vec<T2, cn>() const {
			Vec<T2, cn> v2;
			for (int i = 0; i < channels; ++i) v2[i] = data[i];
			return v2;
		}

		//! dot product computed with the default precision

		inline T dot(const vec& v) const {
			T sum = 0;
			for (int i = 0; i < channels; ++i) sum += v[i] * data[i];
			return sum;
		}

		//! multiply two vectors element-wise

		inline vec mul(const vec& a) const {
			vec v;
			for (int i = 0; i < channels; ++i) v[i] = a[i] * data[i];
			return v;
		}

		//! binary operators

		inline vec operator +(const vec& a) const {
			vec v;
			for (int i = 0; i < channels; ++i) v[i] = data[i] + a[i];
			return v;
		}

		inline vec& operator +=(const vec& a) {
			for (int i = 0; i < channels; ++i) data[i] += a[i];
			return *this;
		}

		inline vec operator -(const vec& a) const {
			vec v;
			for (int i = 0; i < channels; ++i) v[i] = data[i] - a[i];
			return v;
		}

		inline vec& operator -=(const vec& a) {
			for (int i = 0; i < channels; ++i) data[i] -= a[i];
			return *this;
		}

		inline vec operator *(T x) const {
			vec v;
			for (int i = 0; i < channels; ++i) v[i] = data[i] * x;
			return v;
		}

		inline vec& operator *=(T x) {
			for (int i = 0; i < channels; ++i) data[i] *= x;
			return *this;
		}

		//! unary operator

		inline vec operator -() const {
			vec v;
			for (int i = 0; i < channels; ++i) v[i] = -data[i];
			return v;
		}

		//! element access
		inline const T& operator [](int i) const {
			return data[i];
		}
		inline T& operator [](int i) {
			return data[i];
		}

		T data[channels];
	};
	
	typedef Vec<int, 2> Vec2i;
	typedef Vec<int, 3> Vec3i;
	typedef Vec<float, 2> Vec2f;
	typedef Vec<float, 3> Vec3f;
	
	/**
	 * Simple 2d point
	 */
	template <typename T>
	struct Point {
		typedef T scalar;
		typedef Vec<T, 2> vec;
		typedef Point<T> point;
		enum {
			typeDepth = DataDepth<scalar>::value,
			channels = 2,
			type = IM_MAKETYPE(typeDepth, channels),
			exact = DataDepth<scalar>::exact
		};
		Point() : x(0), y(0) {}
		Point(T a, T b) : x(a), y(b) {}
		template <typename T2>
		explicit Point(const Point<T2> &p) : x(p.x), y(p.y) {}
		
		inline T dot(const point &p) const {
			return p.x * x + p.y * y;
		}
		inline T sqNorm() const {
			return dot(*this);
		}
		inline point operator +(const point &p) const {
			return point(x + p.x, y + p.y);
		}
		inline point operator -(const point &p) const {
			return point(x - p.x, y - p.y);
		}
		inline point operator *(T f) const {
			return point(x * f, y * f);
		}
		inline point operator -() const {
			return point(-x, -y);
		}
		inline point mult(const point &p) const {
			return point(x * p.x, y * p.y);
		}
		template <typename I>
		inline point within(const I *image) const {
			return point(
					std::max(T(0), std::min(T(image->width - 1), x)),
					std::max(T(0), std::min(T(image->height - 1), y))
					);
		}
		inline operator Point<int>() const {
			return Point<int>(round(x), round(y));
		}
		inline Point<int> floor() const {
			return Point<int>(std::floor(x), std::floor(y));
		}
		inline Point<int> ceil() const {
			return Point<int>(std::ceil(x), std::ceil(y));
		}
		inline point abs() const {
			return point(std::abs(x), std::abs(y));
		}
		inline bool isOrigin() const {
			return x == 0 && y == 0;
		}
		inline static Point<T> min(const Point<T> &p1, const Point<T> &p2) {
			return Point<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
		}
		inline static Point<T> max(const Point<T> &p1, const Point<T> &p2) {
			return Point<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
		}
		inline T area() const {
			return x * y;
		}
		
		T x;
		T y;
	};
	
	typedef Point<int> Point2i;
	typedef Point<float> Point2f;
	typedef Point<double> Point2d;
    
    template <typename T>
    class Grid {
    public:
        union {
            int rows;
            int height;
        };
        union {
            int cols;
            int width;
        };
        
        Grid() : height(0), width(0), data(){}
        Grid(int h, int w, bool init = false) : height(h), width(w), data() {
            T *ptr = NULL;
            if(init){
                ptr = new T[h * w]();
            } else {
                ptr = new T[h * w];
            }
            if(ptr != NULL) data.reset(ptr);
        }
        
        inline bool empty() const {
            return data;
        }
        
        inline T *ptr() {
            return data.get();
        }
        inline const T *ptr() const {
            
        }
        
        inline T &at(int y, int x) {
            return ptr()[y * width + x];
        }
        inline const T &at(int y, int x) const {
            return ptr()[y * width + x];
        }
        
    private:
        boost::shared_ptr<T> data;
    };
	
	/**
	 * Matrix data pointer
	 */
	typedef boost::shared_ptr<byte> DataPtr;
	
#ifndef SAFE_MAT
#define SAFE_MAT 0
#endif

	/**
	 * Image matrix representation
	 */
	struct Mat {
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
		
		Mat() : flags(IM_UNKNOWN){
		}
		
		Mat(int h, int w, int dataType) : height(h), width(w), flags(dataType){
			int elemSize = IM_SIZEOF(dataType);
			int byteCount = elemSize * h * w;
			if(byteCount > 0){
				byte *content = new byte[byteCount];
				data.reset(content);
				step[0] = elemSize;
				step[1] = w * elemSize;
			} else {
				std::cerr << "No byte for dataType=" << dataType << " with size=" << elemSize << "\n";
				step[0] = step[1] = 0;
			}
		}
		
		inline bool empty() const {
			return !data;
		}
		
		inline static Mat zeros(int rows, int cols, int type) {
			Mat m(rows, cols, type);
			int elemSize = IM_SIZEOF(type);
			int byteCount = elemSize * rows * cols;
			byte *ptr = m.ptr();
			std::fill(ptr, ptr + byteCount, 0);
			return m;
		}
		
		//! Direct pointer access
		inline byte *ptr() {
			return data.get();
		}
		
		//! Pointer access
		template <typename T>
		inline const T *ptr(int y, int x) const {
#if SAFE_MAT
			if(x < 0 || x >= width || y < 0 || y >= height) {
				std::cout << y << "/" << x << "\n";
				mexErrMsgIdAndTxt("MATLAB:img:ptr", "Out of image bounds!");
			}
#endif
			const byte *ref = data.get();
			return reinterpret_cast<const T*>(ref + y * step[1] + x * step[0]);
		}
		template <typename T>
		inline T *ptr(int y, int x) {
#if SAFE_MAT
			if(x < 0 || x >= width || y < 0 || y >= height) {
				std::cout << y << "/" << x << "\n";
				mexErrMsgIdAndTxt("MATLAB:img:ptr", "Ref out of image bounds!");
			}
#endif
			byte *ref = data.get();
			return reinterpret_cast<T*>(ref + y * step[1] + x * step[0]);
		}
		
		//! Element access
		template <typename T>
		inline const T &at(int y, int x) const {
			return *ptr<T>(y, x);
		}
		template <typename T>
		inline T &at(int y, int x) {
			return *ptr<T>(y, x);
		}
		
	private:
		int flags;
		DataPtr data;
		int step[2];
	};
	
	/**
	 * Image type
	 */
	typedef Mat Image;
	typedef Image Mask;

}

#endif	/* ALGEBRA_H */

