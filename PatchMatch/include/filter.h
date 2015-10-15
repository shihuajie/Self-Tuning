/*******************************************************************************
 * filter.h - 2d filter implementations
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef FILTER_H
#define	FILTER_H

#include <cmath>
#include <vector>
#include <iostream>

namespace pm {
	/**
	 * Simple static 2d filter
	 */
	struct Filter {
		Filter() : width(0), weight() {}
		/**
		 * \brief Create a uniform filter
		 * 
         * \param n
		 *			the filter width (size = width * width)
         * @param v
		 *			the uniform weight value
         */
		explicit Filter(int n, float v = 1.0f) : width(n){
			weight.resize(n * n, v);
		}
		// the size
		int width;
		// the weights
		std::vector<float> weight;
		
		// access
		float *operator[](int n){
			return &weight[0] + n * width;
		}
		const float *operator[](int n) const {
			return &weight[0] + n * width;
		}
		
		// normalization
		void normalize(float norm = 1.0f) {
			float sum = 0.0f;
			for (std::vector<float>::iterator it = weight.begin(); it != weight.end(); ++it) {
				sum += *it;
			}
			if(sum < 1e-8) {
				std::cerr << "Warning: using really low filter values, this is dangerous!\n";
				std::cerr << "Consider using a flat filter instead!\n";
			}
			for (std::vector<float>::iterator it = weight.begin(); it != weight.end(); ++it) {
				*it *= norm / sum;
			}
		}
	};
	
	inline Filter flatFilter(int n){
		return Filter(n);
	}
	
	inline Filter gaussianFilter(int N, float sigma) {
		Filter f(N, 0.0f);
		float mid = N * 0.5f;
		float expFactor = -0.5 / sigma;
		for(int y = 0; y < N; ++y){
			float dy = mid - y;
			for(int x = 0; x < N; ++x){
				float dx = mid - x;
				float w = std::exp((dx * dx + dy * dy) * expFactor);
				f[y][x] = w;
			}
		}
		return f;
	}
}

#endif	/* FILTER_H */

