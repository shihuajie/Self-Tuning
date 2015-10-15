/*******************************************************************************
 * ggdt.h - generalized geodesic distance transform
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef GGDT_H
#define	GGDT_H

#include "algebra.h"
#include <cmath>
#include <float.h>

namespace pm {

	template <typename Scalar>
	inline Scalar square(Scalar x){
		return x * x;
	}

	template <typename Scalar, typename Image>
	Image ggdt(const Image &M, Scalar mu, const Image *I = NULL, Scalar delta = 0, int N = 5) {
		// std::cout << "ggdt();\n";
		
		Image D(M.height, M.width, M.type());
		// std::cout << "D(" << M.height << ", " << M.width << ", " << M.type() << ")\n";

		// 1: D_0 = mu M
		for(int y = 0; y < D.height; ++y) {
			for(int x = 0; x < D.width; ++x) {
				D.template at<Scalar>(y, x) = M.template at<Scalar>(y, x) * mu;
			}
		}

		// 2: update using a_k
		//		D(x) = min{ D(x+a_k) + sqrt(|a_k|^2 + delta^2 | I(x) - I(x + a_k) |^2), k=0..4 }
		// const data
		const unsigned int K = 5;
		const int a[K][2] = {
			{ 0, 0 }, { -1, -1 }, { -1, 0 }, { -1, +1}, { 0, -1 }
		};
		const Scalar A[K] = { 0, 2, 1, 2, 1 }; // |a_k|^2
		const Scalar sqrtA[K] = { 0, sqrt(2.0), 1, sqrt(2.0), 1 }; // |a_k|
		const Scalar deltaSq = delta * delta;
		// the iterations
		bool reversed = false;
		N *= 2; // as we need up-left to bottom-right and reversed
		// std::cout << "Looping " << N << " times.\n";
		do {
			// traversal data
			int dir, startX, endX, startY, endY;
			if(!reversed){
				dir = 1;
				startX = 0;
				endX = D.width;
				startY = 0;
				endY = D.height;
			} else {
				dir = -1;
				startX = D.width - 1;
				endX = -1;
				startY = D.height - 1;
				endY = -1;
			}
			
			// traverse the pixels
			if(delta > 0 && I != NULL) {
				for(int y = startY; y != endY; y += dir) {
					for(int x = startX; x != endX; x += dir) {
						// get the minimum of the a_k's
						Scalar minD = D.template at<Scalar>(y, x); // a_0
						for(unsigned int k = 1; k < K; ++k) {
							int y_k = y + dir * a[k][1];
							int x_k = x + dir * a[k][0];
							if(x_k < 0 || y_k < 0 || x_k >= D.width || y_k >= D.height) {
								continue; // we cannot look D[y_k][x_k] as it doesn't exist
							}
							minD = std::min(minD, 
									D.template at<Scalar>(y_k, x_k)
									+ std::sqrt( A[k] + deltaSq
										* square(I->template at<Scalar>(y, x) - I->template at<Scalar>(y_k, x_k) )
									)
							);
						}
						D.template at<Scalar>(y, x) = minD;
					}
				}
			} else {
				for(int y = startY; y != endY; y += dir) {
					for(int x = startX; x != endX; x += dir) {
						// get the minimum of the a_k's
						Scalar minD = D.template at<Scalar>(y, x); // a_0
						for(unsigned int k = 1; k < K; ++k) {
							int y_k = y + dir * a[k][1];
							int x_k = x + dir * a[k][0];
							if(x_k < 0 || y_k < 0 || x_k >= D.width || y_k >= D.height) {
								continue; // we cannot look D[y_k][x_k] as it doesn't exist
							}
							minD = std::min(minD, D.template at<Scalar>(y_k, x_k) + sqrtA[k]);
						}
						D.template at<Scalar>(y, x) = minD;
					}
				}
			}

			// switch the direction
			reversed = !reversed;
		} while(--N >= 0);

		// done
		return D;
	}

}

#endif	/* GGDT_H */

