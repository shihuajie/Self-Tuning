/*******************************************************************************
 * distance.h - abstract patch distance
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef DISTANCE_H
#define DISTANCE_H

#include <limits>
#include <cmath>
#include <patch.h>
#include <gb.h>

namespace pm {

	/**
	 * \brief Distance type
	 */
	enum DistanceType {
		SSD,
		GBASSD,
        CGBASSD,
		wSSD,
        LP
		// ParSSD
	};

	/**
	 * \brief Return whether a distance type uses GBA
	 */
	bool usesGBA(DistanceType type) {
		return type == GBASSD;
	}

	namespace dist {
		
		static bool debug = false;
        
        // #########################################################################
        // ##### SSD ###############################################################
        // #########################################################################

		/**
		 * \brief Simple sum of squared differences
		 */
		template <int channels, typename SourcePatch, typename TargetPatch, typename Scalar>
		Scalar SumSquaredDiff(const Texture *source, const Texture *target,
				const SourcePatch &p1, const TargetPatch &p2, Scalar best) {
			typedef Vec<Scalar, 3> DataType; //modified by huajie 2015-9-6
			if (best < 0) best = std::numeric_limits<Scalar>::max();
			const int width = SourcePatch::width();
			const Scalar invArea = 1.0 / (width * width);
			Scalar sum = 0;
			// assert(SourcePatch::width() == TargetPatch::width());
			for (typename SourcePatch::IndexIterator it = p1.begin(); it; ++it) {
				typename SourcePatch::Index i = *it;
				DataType diff = source->at<DataType>(p1 * i) - target->at<DataType>(p2 * i);
				Scalar d = diff.dot(diff); // squaredNorm<Scalar, channels>(diff);
				sum += d * invArea;
				// sum *  > best || 
				if (!_finite(sum)) return sum;
			}
			return sum;
		}
        
        // #########################################################################
        // ##### L-p Norm ##########################################################
        // #########################################################################

        static const float &pSpace(const float &p = 0.0f){
            static float pValue = 1.0f;
            if(p > 0.0f){
                pValue = p;
            }
            return pValue;
        }
        static const bool &pSpaceRoot(bool doRoot = true) {
            static bool root = true;
            if(!doRoot) root = false;
            return root;
        }
        
		/**
		 * \brief L-p distance
		 */
		template <int channels, typename SourcePatch, typename TargetPatch, typename Scalar>
		Scalar LPDistance(const Texture *source, const Texture *target,
				const SourcePatch &p1, const TargetPatch &p2, Scalar best) {
			typedef Vec<Scalar, channels> DataType;
			if (best < 0) best = std::numeric_limits<Scalar>::max();
			const int width = SourcePatch::width();
			const Scalar invArea = 1.0 / (width * width);
			Scalar sum = 0;
            const float &p = pSpace();
            const bool useRoot = pSpaceRoot();
			// assert(SourcePatch::width() == TargetPatch::width());
			for (typename SourcePatch::IndexIterator it = p1.begin(); it; ++it) {
				typename SourcePatch::Index i = *it;
				DataType diff = source->at<DataType>(p1 * i) - target->at<DataType>(p2 * i);
                Scalar d = 0.0;
                for(int c = 0; c < channels; ++c) {
                    d += std::pow(std::abs(diff[c]), p);
                }
				sum += d * invArea;
				// sum *  > best || 
				if (!_finite(sum)) return sum;
			}
			return useRoot ? std::pow(sum, 1.0f / p) : sum;
		}

		// #########################################################################
		// ##### SSD with Gain and Bias ############################################
		// #########################################################################

		template <int channels, typename SourcePatch, typename TargetPatch, typename Scalar>
		Scalar GainBiasAdjustedSSD(const Texture *source, const Texture *target,
				const SourcePatch &p1, const TargetPatch &p2, Scalar best) {
			typedef Vec<Scalar, channels> DataType;
			typedef Vec<Scalar, 3> Vec3;
			if (best < 0) best = std::numeric_limits<Scalar>::max();
			const Scalar invArea = 1.0 / (SourcePatch::width() * SourcePatch::width());

			// bias and gain
			Vec3 gain, bias;
			GainBias<Scalar>::template compute<channels>(source, target, p1, p2, gain, bias);

			// the distance sum
			Scalar sum = 0;
			// assert(SourcePatch::width() == TargetPatch::width());
			for (typename SourcePatch::IndexIterator it = p1.begin(); it; ++it) {
				typename SourcePatch::Index i = *it;
				DataType q = source->at<DataType>(p1 * i);
				DataType p = target->at<DataType>(p2 * i);

				// gain/bias adjustment of p
				GainBias<Scalar>::template applyOn<channels>(p, gain, bias);

				// different and contribution to the distance sum
				DataType diff = q - p;
				Scalar d = diff.dot(diff) * invArea;
				sum += d;
				// sum *  > best || 
				if (!_finite(sum)) return sum;
			}
			return sum;
		}
        
        template <int channels, typename SourcePatch, typename TargetPatch, typename Scalar>
		Scalar CachedGBA_SSD(const Texture *source, const Texture *target,
				const SourcePatch &p1, const TargetPatch &p2, Scalar best) {
			typedef Vec<Scalar, channels> DataType;
			typedef Vec<Scalar, 3> Vec3;
			if (best < 0) best = std::numeric_limits<Scalar>::max();
			const Scalar invArea = 1.0 / (SourcePatch::width() * SourcePatch::width());

			// bias and gain
			Vec3 gain, bias;
			GainBias<Scalar>::fetch(GainBias<Scalar>::Target, GainBias<Scalar>::Exemplar, p1, p2, gain, bias);

			// the distance sum
			Scalar sum = 0;
			// assert(SourcePatch::width() == TargetPatch::width());
			for (typename SourcePatch::IndexIterator it = p1.begin(); it; ++it) {
				typename SourcePatch::Index i = *it;
				DataType q = source->at<DataType>(p1 * i);
				DataType p = target->at<DataType>(p2 * i);

				// gain/bias adjustment of p
				GainBias<Scalar>::template applyOn<channels>(p, gain, bias);

				// different and contribution to the distance sum
				DataType diff = q - p;
				Scalar d = diff.dot(diff) * invArea;
				sum += d;
				// sum *  > best || 
				if (!_finite(sum)) return sum;
			}
			return sum;
		}
		
		// #########################################################################
		// ##### Weighted SSD ######################################################
		// #########################################################################
		
		Texture &ssdWeight(Texture *w = NULL) {
			static Texture ssd_weight;
			if(w != NULL) {
				ssd_weight = *w;
			}
			return ssd_weight;
		}
		
		/**
		 * \brief Weighted sum of squared differences
		 */
		template <int channels, typename SourcePatch, typename TargetPatch, typename Scalar>
		Scalar WeightedSumSquaredDiff(const Texture *source, const Texture *target,
				const SourcePatch &p1, const TargetPatch &p2, Scalar best) {
			typedef Vec<Scalar, channels> DataType;
			if (best < 0) best = std::numeric_limits<Scalar>::max();
			const int width = SourcePatch::width();
			const Scalar invArea = 1.0 / (width * width);
			Scalar sum = 0;
			// weight
			Texture &weight = ssdWeight();
			Scalar wSum = 0.0;
			// assert(SourcePatch::width() == TargetPatch::width());
			for (typename SourcePatch::IndexIterator it = p1.begin(); it; ++it) {
				typename SourcePatch::Index i = *it;
				typename SourcePatch::PixLoc l1 = p1 * i;
				DataType diff = source->at<DataType>(p1 * i) - target->at<DataType>(p2 * i);
				Scalar w = weight.at<Scalar>(l1); // the weight in the source
				Scalar d = diff.dot(diff) * w * w;
				wSum += w * w;
				sum += d;
				// sum *  > best || 
				if (!_finite(sum)) return sum;
			}
			return sum / (wSum > 0 ? wSum : 1.0);
		}

	}

	/**
	 * \brief Distance factory and type helper
	 */
	template <typename Patch, typename Scalar>
	struct Distance {
		/// Source patch type
		typedef typename Patch::OriginalPatchType SourcePatch;
		/// Target patch type
		typedef Patch TargetPatch;
		/**
		 * Dist function pointer type
		 */
		typedef Scalar(*Function)(const Texture *, const Texture *, const SourcePatch &, const TargetPatch &, Scalar);

		/**
		 * \brief Return the distance corresponding to a given type
		 * 
		 * \param type the distance type looked for
		 * \return the corresponding distance function
		 */
		template <int Channels>
				static Function get(DistanceType type) {
			switch (type) {
				case wSSD:
					if(!dist::ssdWeight().empty()) {
						return &dist::WeightedSumSquaredDiff<Channels, SourcePatch, TargetPatch, Scalar>;
					}
					std::cerr << "Empty ssd weight!\n";
				case SSD: return &dist::SumSquaredDiff<Channels, SourcePatch, TargetPatch, Scalar>;
                case LP: return &dist::LPDistance<Channels, SourcePatch, TargetPatch, Scalar>;
				case GBASSD: return &dist::GainBiasAdjustedSSD<Channels, SourcePatch, TargetPatch, Scalar>;
                case CGBASSD:
                    if(GainBias<Scalar>::cached()){
                        std::cerr << "Caches have not been loaded!\n";
                        return NULL;
                    }
                    return &dist::GainBiasAdjustedSSD<Channels, SourcePatch, TargetPatch, Scalar>;
				default:
					std::cerr << "Invalid distance type!";
					return NULL;
			}
		}
	};

}

#endif