/*******************************************************************************
 * histogram.h - histogram algorithms for voting
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef HISTOGRAM_H
#define	HISTOGRAM_H

#include <patch.h>
#include <vector>

namespace pm {
	
	/**
	 * Channel range
	 */
	struct Range {
		float from, to;
		Range(float a, float b) : from(a), to(b) {}
		Range() : from(0.0f), to(0.0f) {}
		
		inline float size() const {
			return to - from;
		}
	};
	
	namespace voting {
		
		template <typename Patch, typename Scalar = float>
		class Histogram {
		public:
			Histogram() {}
			Histogram(const Image *image, int channel, int numBins, const Range &r)
			: bins(numBins, 0), range(r), binCount(numBins), binSize(r.size() / numBins){
				// populate the histogram from an image channel
				for(int y = 0; y < image->rows; ++y) {
					for(int x = 0; x < image->cols; ++x) {
						const Scalar *ptr = image->ptr<Scalar>(y, x);
						Scalar v = ptr[channel];
						int bi = binOf(v);
						bins[bi] += 1;
					}
				}
			}
			
			/// Return the bin index corresponding to a given value
			inline int binOf(float value) const {
				int idx = std::floor(float(value - range.from) / binSize);
				return std::max(0, std::min(idx, binCount - 1));
			}
			
			/// Access operators
			
			inline int count(Scalar v) const {
				return bins[binOf(v)];
			}
			
			inline int &count(Scalar v) {
				return bins[binOf(v)];
			}
			
		private:
			std::vector<int> bins;
			Range range;
			int binCount;
			float binSize;
		};
	}
}

#endif	/* HISTOGRAM_H */

