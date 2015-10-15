/*******************************************************************************
 * rng.h - random number utils
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef RNG_H
#define	RNG_H

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

namespace pm {
	// Random Number Generator
	typedef float (*RNG)(void);

	float unif01(void) {
		return float(rand()) / float(RAND_MAX);
	}
	
	/**
	 * \brief Initialize the RNG with a given seed
	 * 
     * \param s the seed to use
     */
	inline void seed(unsigned int s) {
		srand(s);
	}
	
	/**
	 * \brief Initialize the RNG with the current time (seconds!)
     */
	inline void seed() {
		seed(time(NULL));
	}

	/// Discrete Uniform RV ~ Unif[a, b]

	inline int uniform(RNG rand, int a, int b) {
		if (a >= b) std::cout << "Invalid uniform: (" << a << " to " << b << ")\n";
		int r = a + std::floor(rand() * (b - a + 1));
		if (r > b) return b;
		return r;
	}

	inline float uniform(RNG rand, float a, float b) {
		if (a >= b) std::cout << "Invalid uniform: [" << a << " to " << b << "]\n";
		return a + (b - a) * rand();
	}

	/// Discrete Bernoulli RV ~ Bernoulli(p)

	inline bool bernoulli(RNG rand, float p = 0.5f) {
		return rand() <= p;
	}
	
	/// Gaussian Noise
	template <typename T>
	inline T gaussian(RNG rand, T sigma, bool noStore = false) {
		// Box-Muller transform
		// @see http://projecteuclid.org/DPubS?verb=Display&version=1.0&service=UI&handle=euclid.aoms/1177706645&page=record
		// @see http://en.wikipedia.org/wiki/Box-Muller_transform
		
		// the last instances
		static T u1; // radius
		static T u2; // angle
		static bool hasSpare = false; // whether we have only used them once

		// reuse if they're fresh
		if(hasSpare){
			hasSpare = false;
			return sigma * std::sqrt(u1) * std::sin(u2);
		}

		// renew data
		u1 = rand();
		if(u1 < 1e-100) u1 = 1e-100;
		u1 = -2 * std::log(u1); // squared radius, inv-exponentially distributed
		u2 = rand() * M_PI * 2; // angle, uniformly distributed over [0;2pi]

		// available for afterwards
		hasSpare = !noStore; //< unless we asked not to keep it
		// gaussian version
		return sigma * std::sqrt(u1) * std::cos(u2);
	}
	
	template <typename Point, typename T>
	inline Point gaussian2d(RNG rand, T sigma) {
		// squared radius, inv-exponentially distributed
		T u1 = rand();
		if(u1 < 1e-100) u1 = 1e-100;
		u1 = -2 * std::log(u1);
		
		// angle, uniformly distributed over [0;2pi]
		T u2 = rand() * M_PI * 2;
		
		// our sample
		T r_s = std::sqrt(u1) * sigma;
		return Point(r_s * std::cos(u2), r_s * std::sin(u2)); 
	}

	/**
	 * \brief Knuth in-place shuffling
	 * 
	 * \param r
	 *			the random number source
	 * \param index
	 *			the index to shuffle
	 * \param N
	 *			the size of the index
	 * \see http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
	 * \note randomness is limited as the PRNG is now sophisticated, 
	 * thus the shuffle will be bad and biased, but it should be enough ...
	 */
	template <typename T>
	inline void knuth_shuffle(RNG r, T *index, int N) {
		for (int i = N - 1; i > 0; --i) {
			int j = uniform(r, 0, i);
			std::swap(index[j], index[i]);
		}
	}
}

#endif	/* RNG_H */

