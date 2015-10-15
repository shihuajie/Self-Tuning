/*******************************************************************************
 * patch.h - abstract image patch
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef PATCH_H
#define PATCH_H

#include <cmath>
#include <cstdlib>

#ifndef M_PI
#define M_PI            3.14159265358979323846
#define M_PI_4          0.785398163397448309616
#endif

#include <rng.h>
#include <texture.h>

#define DEFAULT_PATCH_SIZE 7

namespace pm {

	// Patch type

	enum PatchType {
		Basic,
		BasicFloat,
		Affine,
		BasicMultiRes,
		BasicFloatMultiRes,
		AffineMultiRes
		// BilinearAffine
	};
	
	/**
	 * \brief Grid index iterator
	 */
	template <typename Index>
	struct LinearPatchIterator {
		typedef LinearPatchIterator<Index> IndexIterator;
		Index i;
		const int width;

		LinearPatchIterator(int y, int x, int w) : i(x, y), width(w){}

		inline IndexIterator &operator ++(){
			++i.x;
			if(i.x >= width){
				++i.y;
				i.x = 0;
			}
			return *this;
		}

		inline const Index &operator *() const {
			return i;
		}
		inline operator bool() const {
			return i.y < width;
		}
	};
	
	/**
	 * Patch base type for iterations
     */
	template <typename Patch>
	struct BasicGrid {
		
		static int depth(int = 0) {
			mexErrMsgIdAndTxt("MATLAB:basicgrid:depth", "Called depth on a simple resolution patch!");
			return 0;
		}
		
		/// The type of the index iterator
		typedef Point<int> Index;
		typedef LinearPatchIterator<Index> IndexIterator;
		
		IndexIterator begin() const {
			return IndexIterator(0, 0, Patch::width());
		}
		const IndexIterator end() const {
			int w = Patch::width();
			return IndexIterator(w, 0, w);
		}
	};

	namespace patch {

		/**
		 * \brief Create a random patch
		 * 
		 * \param rand
		 *			the random number generator
		 * \param parent
		 *			the image parent of the patch
		 * \param patch
		 *			the patch to randomly initialize
		 */
		template <typename Patch>
		void randomInit(RNG rand, const Image *parent, Patch &patch);
		
		/**
		 * \brief Initialize a new random patch taken from an old one
		 * 
		 * \param rand
		 *			the random number generator
		 * \param parent
		 *			the image parent of the patches
		 * \param oldPatch
		 *			the patch we randomly generate from
		 * \param newPatch
		 *			the new patch we randomly initialize
		 * \param windowSize
		 *			the random search window size
		 * \return whether the patch is valid (within the frame)
		 */
		template <typename Patch>
		bool random(RNG rand, const Image *parent,
				const Patch &oldPatch, Patch &newPatch, int windowSize);
        
        /**
		 * \brief Initialize a new random patch taken from an old one using aligned search
		 * 
		 * \param rand
		 *			the random number generator
		 * \param parent
		 *			the image parent of the patches
		 * \param oldPatch
		 *			the patch we randomly generate from
		 * \param newPatch
		 *			the new patch we randomly initialize
		 * \param g1
         *          the first direction
         * \param g2
         *          the second direction
		 * \return whether the patch is valid (within the frame)
		 */
		template <typename Patch>
		bool aligned(RNG rand, const Image *parent,
				const Patch &oldPatch, Patch &newPatch, const Point2f &g1, const Point2f &g2, float jitter);

		/**
		 * \brief Initialize a delta patch
		 * 
		 * \param patch
		 *			the fixed patch
		 * \param delta
		 *			the delta patch to initialize
		 * \param dy
		 *			the y delta
		 * \param dx
		 *			the x delta
		 */
		template <typename Patch>
		void deltaPatch(const Patch &patch, Patch &delta, int dy, int dx);

		/**
		 * \brief Check whether a patch is within a given frame
		 * 
		 * \param patch
		 *			the patch
		 * \param maxY
		 *			the frame Y size
		 * \param maxX
		 *			the frame X size
		 * \return whether the patch is within the frame (0, 0, maxY, maxX)
		 */
		template <typename Patch>
		bool isWithin(const Image *frame, const Patch &patch);

		/**
		 * \brief Compute the [0-1] coherence between two patches
		 * 
		 * \param p1
		 *			the first patch
		 * \param p2
		 *			the second patch
		 * \param dy
		 *			the expected dy from p1 to p2
		 * \param dx
		 *			the expected dx from p1 to p2
		 * \return 0 if not coherent, in (0;1] when coherent, 1 if fully coherent
		 */
		template <typename Patch>
		typename Patch::Coherence coherence(const Patch &p1, const Patch &p2, int dy, int dx);
	}
	
	// displacement of location
	typedef int PatchDisplacement;
	
	template <typename Patch>
	inline bool random(RNG rand, const Image *parent,
			const Patch &oldPatch, Patch &newPatch,
			int windowSize, PatchDisplacement minDisp, int y0, int x0) {
		using patch::random;
		if(minDisp <= 0) return random(rand, parent, oldPatch, newPatch, windowSize);
			// draw until we find a correct one otherwise
		int i = 0;
		bool valid = false;
		PatchDisplacement d = 0;
		do {
			// TODO rejection sampling is maybe not the best! 
			valid = random(rand, parent, oldPatch, newPatch, windowSize);
			d = std::abs(newPatch.x - x0) + std::abs(newPatch.y - y0); // L1 distance
		} while(++i < 100 && d < minDisp);
		return valid;
	}
    
    template <typename Patch>
	inline bool aligned(RNG rand, const Image *parent,
			const Patch &oldPatch, Patch &newPatch,
            const Point2f &g1, const Point2f &g2, float jitter,
			PatchDisplacement minDisp, int y0, int x0) {
		using patch::aligned;
		if(minDisp <= 0) return aligned(rand, parent, oldPatch, newPatch, g1, g2, jitter);
			// draw until we find a correct one otherwise
		int i = 0;
		bool valid = false;
		PatchDisplacement d = 0;
		do {
			// TODO rejection sampling is maybe not the best! 
			valid = aligned(rand, parent, oldPatch, newPatch, g1, g2, jitter);
			d = std::abs(newPatch.x - x0) + std::abs(newPatch.y - y0); // L1 distance
		} while(++i < 100 && d < minDisp);
		return valid;
	}
}

#endif