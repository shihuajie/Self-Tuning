/*******************************************************************************
 * multires.h - multi-resolution patches
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef MULTIRES_H
#define	MULTIRES_H

#include <patch.h>
#include <rng.h>
#include <vector>

namespace pm {

	template <typename Patch>
	struct MultiRes : public Patch {
		typedef MultiRes<typename Patch::OriginalPatchType> OriginalPatchType;
		typedef typename Patch::Index ParentIndex;
		typedef MultiDepthPoint<int> Index;
		typedef typename Patch::PixLoc ParentPixLoc;
		typedef MultiDepthPoint<typename ParentPixLoc::scalar> PixLoc;
		
		MultiRes() : Patch() {}
		MultiRes(int y, int x) : Patch(y, x) {}
		
		bool operator==(const MultiRes<Patch> &p) const {
			return Patch::operator ==(reinterpret_cast<const Patch&>(p));
		}

		static int depth(int numLevels = 0) {
			static int num = 0; //< the number of different resolution levels (in addition to the normal one)
			if (numLevels == 0) return num;
			// else we reset it when it's valid
			if (numLevels > 0) {
				num = numLevels;
			}
		}

		struct MultiDepthIterator {

			MultiDepthIterator(int y, int x, int d, int depth, int width) : i(y, x, d), depth(depth), width(width) {
				if(d > 0) {
					const int d0 = d;
					int F = 1;
					while(--d > 0) F *= 3;
					// => F = 3**(d0-1)
					P = F * width;
					TL = Index(1, 1, d0) * ((F * 3 - 1) / 2); 
				} else {
					TL = Index(0, 0, 0);
					P = 1;
				}
			}

			/**
			 * TL(0) = (0,0)
			 * TL(1) = -width (1,1)
			 * TL(n) = TL(n-1) - P(n) = -width * (3**N - 1) / 2
			 */
			Index TL;
			/**
			 * P(0) = 1
			 * P(1) = width
			 * P(n) = 3 * P(n-1) = 3**(n-1) * width
			 */
			int P;
			Index i;
			const int depth, width;

			inline MultiDepthIterator &operator ++() {
				if (i.depth == 0) {
					++i.x;
					if (i.x >= width) {
						++i.y;
						i.x = 0;
					}
					if (i.y >= width) {
						i = Index(0, 0, 1);
						TL = Index(1, 1, 1) * (-width);
						P = width;
					}
				} else {
					++i.x;
					if(i.x > 2) {
						++i.y;
						i.x = 0;
					} else if(i.x == 1 && i.y == 1){
						i.x = 2; // we skip (1, 1, d)
					}
					// next depth level?
					if(i.y > 2) {
						i = Index(0, 0, i.depth + 1);
						P *= 3; // P(n) = 3 * P(n-1)
						TL = Index(TL.x - P, TL.y - P, TL.depth + 1); // TL(n) = TL(n-1) - P(n) * (1,1)
					}
				}
				return *this;
			}

			inline const Index operator *() const {
				Index ind = TL + i * P;
				return ind;
			}
			inline operator bool() const {
				return i.depth <= depth; // equality is ok as we'd be in the last depth level
			}
		};
		
		//! Pixel transformation
		inline PixLoc transform(const Index &i) const {
			return PixLoc(Patch::transform(i.y, i.x), i.depth);
		}
		inline PixLoc operator *(const Index &i) const {
			return transform(i);
		}
		inline ParentPixLoc transform(int py, int px) const {
			return Patch::transform(py, px);
		}
		inline ParentPixLoc transform(const ParentIndex &i) const {
			return Patch::transform(i.y, i.x);
		}
		inline ParentPixLoc operator *(const ParentIndex &i) const {
			return Patch::transform(i);
		}
		
		typedef MultiDepthIterator IndexIterator;
		
		MultiDepthIterator begin() const {
			return MultiDepthIterator(0, 0, 0, depth(), Patch::width());
		}
		const MultiDepthIterator end() const {
			int d = depth();
			return MultiDepthIterator(-1, -1, d + 1, d, Patch::width());
		}
	};
	
	template <typename Patch>
	inline void randomInit(RNG rand, const Image *parent, MultiRes<Patch> &patch) {
		randomInit(rand, parent, reinterpret_cast<Patch&>(patch));
	}

	template <typename Patch>
	inline bool random(RNG rand, const Image *parent, const MultiRes<Patch> &oldPatch,
			MultiRes<Patch> &newPatch, int windowSize) {
		return random(rand, parent, reinterpret_cast<const Patch&>(oldPatch), 
				reinterpret_cast<Patch&>(newPatch), windowSize);
	}
    
    template <typename Patch>
	inline bool aligned(RNG rand, const Image *parent,
			const MultiRes<Patch> &oldPatch,
            MultiRes<Patch> &newPatch, 
            const Point2f &g1, const Point2f &g2, float jitter) {
        return aligned(rand, parent, reinterpret_cast<const Patch&>(oldPatch), 
				reinterpret_cast<Patch&>(newPatch), g1, g2, jitter);
    }

	template <typename Patch>
	inline void deltaPatch(const MultiRes<Patch> &patch, MultiRes<Patch> &delta,
			int dy, int dx) {
		deltaPatch(reinterpret_cast<const Patch&>(patch), reinterpret_cast<Patch&>(delta), dy, dx);
	}

	template <typename Patch>
	inline bool isWithin(const Image *parent, const MultiRes<Patch> &patch) {
		return isWithin(parent, reinterpret_cast<const Patch&>(patch));
	}

	/// Coherence: 1.0f if x'=x+dx and y'=y+dy, 0.0f otherwise
	template <typename Patch>
	inline typename Patch::Coherence coherence(const MultiRes<Patch> &p1, const MultiRes<Patch> &p2, int dy, int dx){
		return coherence(reinterpret_cast<const Patch&>(p1), reinterpret_cast<const Patch&>(p2), dy, dx);
	}
}

#endif	/* MULTIRES_H */

