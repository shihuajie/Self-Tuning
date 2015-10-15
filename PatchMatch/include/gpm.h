/*******************************************************************************
 * gpm.h - generalized patch match algorithm
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef GPM_H
#define GPM_H

#include <boost/shared_ptr.hpp>
#include <distance.h>
#include <map>
#include <nnf.h>
#include <mexutil.h>
#include <segment.h>
#include <set>
#include <vector>

namespace pm {

	struct ConvergenceData {
		int propCount;
		int rsCount;
        int asCount;
		int isCount;
		float maxDist;
		float meanDist;
		float cohRatio;
		float occRatio;
	};

	enum IterationOrder {
		Normal = 0,
		Reverse = 1
	};

	/**
	 * Parameters for the NN computations
	 */
	struct NNSettings {
		// -- patch
		int patchSize;
		int multires;
		// -- iterations
		int iterations;
		int untilConvergence;
		IterationOrder startOrder;
		// -- mask
		Image mask;
		bool scramble;
		// -- comp
		CompletenessParameters completeness;
		bool completePropagation;
		// -- random search
		bool coherentRandSearch;
		int randSearch;
		int maxRandSearch;
		bool mixedRandSearch;
		int windowSize; // int minWindowSize, maxWindowSize;
		bool windowSizeDecr;
        // -- aligned search
        int alignedSearch;
		int maxAlignedSearch;
        Point2f alignG1, alignG2;
        float alignP1, alignP2;
        float alignJitter;
		// -- incomplete search
		int incompleteSearch;
		int maxIncompleteSearch;
		int jumpBufferSize;
		float jumpSigmaRatio;
		// -- distance
		DistanceType distType;
		// -- identity distance
		PatchDisplacement minPatchDisp;
		// -- output
		bool storeConvergence, storeOccRatio;
		std::vector<ConvergenceData> convergence;
        bool occDist;
        // weight dataset sequence
        float weightIndex; //add by huajie 2015-9-11
		NNSettings() : mask(), completeness(), convergence() {
			patchSize = 10;
			multires = 0;
			iterations = 6;
			untilConvergence = 0;
			startOrder = Normal;
			scramble = false;
			completePropagation = false;
			coherentRandSearch = false;
			randSearch = 6;
			maxRandSearch = 6;
			mixedRandSearch = true;
			windowSize = 0;
			windowSizeDecr = true;
            alignedSearch = 0;
			maxAlignedSearch = 6;
            alignP1 = alignP2 = 0.0f;
            alignJitter = 0.0f;
			incompleteSearch = 0;
			maxIncompleteSearch = 6;
			jumpBufferSize = 0;
			jumpSigmaRatio = 0.42f;
			distType = SSD;
			minPatchDisp = 0;
			storeConvergence = false;
            storeOccRatio = false;
            occDist = false;
            weightIndex = 0.0; //add by huajie 2015-9-11
		}
	};

	template <typename Patch, typename Scalar>
	NearestNeighborField<Patch, Scalar> *nnf(
			const Texture *source, 
			const Texture *target,
			NearestNeighborField<Patch, Scalar> *prev = NULL,
			typename NearestNeighborField<Patch, Scalar>::Extension *ext = NULL, 
			NNSettings &settings = NNSettings());
	
	struct PatchSegment{
		int size;
		Point2i offset;
		std::set<unsigned int> neighbors;
		bool valid;
		PatchSegment *parent;
		
		PatchSegment() : size(0), neighbors(), valid(true), parent(NULL) {
		}
		
		bool isRoot() const {
			return parent == NULL;
		}
		
		PatchSegment *root(){
			if(isRoot()) return this;
			return parent = parent->root();
		}
	};

	/**
	 * \brief Compute the NNF for a fixed-N-channel source and target
	 * 
	 * \param source
	 *			the source image
	 * \param target
	 *			the target image
	 * \param prev
	 *			the previous nnf to use
	 * \param extNNF
	 *			the nnf extension to use
	 * \param settings
	 *			the algorithm parameters
	 * \return the computed nnf
	 */
	template <typename Patch, typename Scalar, int channels>
	NearestNeighborField<Patch, Scalar> *nnf_n(const Texture *source, const Texture *target,
			NearestNeighborField<Patch, Scalar> *prev,
			typename NearestNeighborField<Patch, Scalar>::Extension *extNNF, 
			NNSettings &settings) {
		typedef NearestNeighborField<Patch, Scalar> NNF;
		// create the required nnf field
		NNF *field = NULL;
		if (prev) {
			field = prev;
            field->weightIndex = settings.weightIndex;
			field->compParams = settings.completeness;
			field->distFunc = Distance<Patch, Scalar>::template get<channels>(settings.distType);
			field->minPatchDisp = settings.minPatchDisp;
			if (!field->check()) return NULL; // invalid nnfs should not happen!
			field->calcDistances(); //< most likely, we want to recompute the distances now
			field->calcCompleteness();
			// /!\ if not, then ... we're doing again a new run? why?
			if (settings.scramble) field->scramble(settings.mask);
		} else {
			field = new NNF(source, target, true,
					Distance<Patch, Scalar>::template get<channels>(settings.distType),
					settings.completeness, settings.minPatchDisp);
			field->randomize();
		}
		if(extNNF != NULL) {
			field->useExtension(extNNF);
		}
		if(settings.incompleteSearch > 0) {
			field->initJumpBuffer(settings.jumpBufferSize, settings.jumpSigmaRatio);
		}

		std::cout << "NNF: " << field->height << " x " << field->width << "\n";
		// real window size
		if (settings.windowSize <= 0) {
			settings.windowSize = std::max(target->cols, target->rows);
		}
		// convergence data
		if (settings.storeConvergence)
			settings.convergence.resize(settings.iterations);
		// segment data
		typedef Segmentation<NoData> Segments;
		boost::shared_ptr<Segments> segmentPtr;
		// alternate between scanline and reverse-scanline order
		bool reverseOrder = settings.startOrder == Reverse;
		bool doProp = true;
		int numRandSearch = settings.randSearch;
		int numAlignedSearch = settings.alignedSearch;
		int numIncompleteSearch = settings.incompleteSearch;
		int numProp = settings.iterations;
		Image &mask = settings.mask;
		for (int i = 0; i < numProp; i++, reverseOrder = !reverseOrder) {
			int dx, dy, startX, startY, endX, endY;
			if (reverseOrder) {
				dx = dy = -1;
				startX = field->width - 1;
				startY = field->height - 1;
				endX = -1;
				endY = -1;
			} else {
				dx = dy = 1;
				startX = 0;
				startY = 0;
				endX = field->width;
				endY = field->height;
			}

			// abort propagation if there is no need to do it
			if (!doProp && extNNF == NULL && !settings.mixedRandSearch) endY = startY;
			doProp = false; // by default, we don't need to propagate more
			int propCount = 0, simCount = 0, rsCount = 0, asCount = 0, isCount = 0;

			// for each patch, do a propagation and random search step
			for (int y = startY; y != endY; y += dy) {
				for (int x = startX; x != endX; x += dx) {
					// check the mask
					if (!mask.empty() && mask.at<float>(y, x) >= 1.0f) {
						continue; // we shouldn't touch that nnf position
					}
					// check whether propagation makes sense
					bool fullyCoherent = field->coherence(y, x) >= 4.0f;
					if(!fullyCoherent){
						////////////////////////////////////////////////////////////
						// coherence ///////////////////////////////////////////////
						////////////////////////////////////////////////////////////
						// propagate surrounding patches to the current one
						// Note: surrounding means the previous ones, else we're
						// destroying our work at each new step of the propagation
						bool didProp = field->propagation(y, x, dy, dx, settings.completePropagation);

						// conv data
						if (didProp) {
							doProp = true;
							++propCount; // one more instance of valid propagation
						}
					}
					if (!fullyCoherent || settings.coherentRandSearch){
					
						////////////////////////////////////////////////////////////
						// similarity //////////////////////////////////////////////
						////////////////////////////////////////////////////////////
						if(extNNF != NULL) {
							bool didSim = field->similarityPropagation(y, x, dy, dx, settings.completePropagation);

							// conv data
							if(didSim) {
								doProp = true;
								++simCount;
							}

							// switch to the next one
							// XXX should we switch even if the last qx succeeded?
							field->nextSimilarity(y, x);
						}

						////////////////////////////////////////////////////////////
						// random search ///////////////////////////////////////////
						////////////////////////////////////////////////////////////

						// mixed random search
						if (settings.mixedRandSearch && (!fullyCoherent || settings.coherentRandSearch)) {
							bool res = false;
							for (int r = 0; r < numRandSearch; ++r) {
								// do a random search to get better solutions
								res = field->randomSearch(y, x, settings.windowSize, settings.maxRandSearch);
							}
							if (res) {
								doProp = true;
								++rsCount;
							}

							// aligned search?
							res = false;
							for(int r = 0; r < numAlignedSearch; ++r) {
								// do an aligned search
								res = field->alignedSearch(y, x,
										settings.alignG1, settings.alignG2, settings.alignJitter, 
										settings.maxAlignedSearch);
							}
							if (res) {
								doProp = true;
								++asCount;
							}

							// incomplete search?
							res = false;
							for(int r = 0; r < numIncompleteSearch; ++r) {
								// do an aligned search
								res = field->incompleteSearch(y, x, settings.maxIncompleteSearch);
							}
							if (res) {
								doProp = true;
								++isCount;
							}
						}
					}
				}
			}

			// separate random search
			if(!settings.mixedRandSearch && (numRandSearch > 0 || numAlignedSearch > 0 || numIncompleteSearch > 0)){
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
				for (int x = 0; x < field->width; ++x) {
					for (int y = 0; y < field->height; ++y) {
						if (!settings.coherentRandSearch && field->coherence(field->get(y, x), y, x) >= 4) {
							continue; // we do not search for fully coherent patches
						}
						// random search
						bool res = false;
						for (int r = 0; r < numRandSearch; ++r) {
							// do a random search to get better solutions
							res = field->randomSearch(y, x, settings.windowSize, settings.maxRandSearch);
						}
						if (res) {
							doProp = true;
#if _OPENMP
#pragma	omp atomic
#endif
							++rsCount;
						}
						// aligned search
						res = false;
						for (int r = 0; r < numAlignedSearch; ++r) {
							// do a random search to get better solutions
							res = field->alignedSearch(y, x,
									settings.alignG1, settings.alignG2, settings.alignJitter, 
									settings.maxAlignedSearch);
						}
						if (res) {
							doProp = true;
#if _OPENMP
#pragma	omp atomic
#endif
							++asCount;
						}
						// incomplete search
						res = false;
						for(int r = 0; r < numIncompleteSearch; ++r) {
							// do an aligned search
							res = field->incompleteSearch(y, x, settings.maxIncompleteSearch);
						}
						if (res) {
							doProp = true;
#if _OPENMP
#pragma	omp atomic
#endif
							++isCount;
						}
					}
				}
			}
			// reduction of the window size
			if (settings.windowSizeDecr && settings.windowSize > 5) {
				settings.windowSize >>= 1; // divide by 2!
			}

			// last iteration?
			if (i == numProp - 1) {
				if (doProp) {
					numProp = std::min(numProp + 1, settings.untilConvergence); // one more propagation
					numRandSearch = 0; // no more random search
					numAlignedSearch = 0;
					numIncompleteSearch = 0;
				}
				// when storing convergence, we have to extend the container
				if (settings.storeConvergence && i >= settings.iterations) {
					settings.convergence.push_back(ConvergenceData());
				}
			} else {
				// if nothing changed, we may want to search more
				if (!doProp) {
					++numRandSearch;
				} else {
					// else we can search once less, but at least once
					numRandSearch = std::max(1, numRandSearch - 1);
					// XXX or do we keep searching more?
					// XXX or do we go back to numRandSearch=1 ?
					// XXX or do we half numRandSearch?
				}
			}
			// should we invalidate the number of random searches to be done?
			if (settings.randSearch <= 0) numRandSearch = 0;

			// for debug purposes
			Scalar maxD = field->maxDistance();
			Scalar meanD = field->meanDistance();
			Scalar cohRatio = field->coherence() / field->size;
			Scalar occRatio = settings.storeOccRatio ? field->meanPatchOccRatio() : 0;
			if (settings.storeConvergence) {
				ConvergenceData &cv = settings.convergence[i];
				cv.propCount = propCount;
				cv.rsCount = rsCount;
				cv.asCount = asCount;
				cv.isCount = isCount;
				cv.maxDist = maxD;
				cv.meanDist = meanD;
				cv.cohRatio = cohRatio;
				cv.occRatio = occRatio;
			}
			std::cout << (i + 1) << ". iteration completed [p=" << propCount;
			std::cout << ", rs=" << rsCount << ", sc=" << simCount;
			std::cout << ", as=" << asCount << ", ic=" << isCount << "]";
			if(settings.incompleteSearch > 0){
				int vsc, tsc, sjc, src;
				field->getSampleCount(vsc, tsc, sjc, src);
				std::cout << "{samp " << int(vsc * 100.0f / tsc) << "%, jump ";
				std::cout << int(sjc * 100.0f / vsc) << "%, reset " << src << "}";
			}
			std::cout << "(maxDist=" << maxD << ", meanDist=" << meanD;
			std::cout << ", coh=" << cohRatio << ", occRatio=" << occRatio << ").\n";
			if (!_finite(maxD) || !_finite(meanD)) {
				for (int y = 0; y < field->height; ++y) {
					for (int x = 0; x < field->width; ++x) {
						std::cout << "@" << y << "/" << x << ": ";
						std::cout << field->distance(y, x) << " =?= ";
						std::cout << field->distance(y, x, field->get(y, x)) << "\n";
					}
				}
				mexErrMsgIdAndTxt("MATLAB:gpm:invalid_distances", "The distances are wrong!");
			}
			
			// is it useful to go farther?
			if(!doProp && settings.randSearch <= 0) {
				return field;
			}
			
			// post-iteration operations
			if(numIncompleteSearch > 0) {
				field->updateJumpBuffer();
			}
		}

		return field;
	}

	// #########################################################################
	// ##### Multi-channel implementation ######################################
	// #########################################################################

	template <typename Patch, typename Scalar, int channels = 1 >
	struct MultiChannelNNF {

		inline NearestNeighborField<Patch, Scalar> *operator()(
				const Texture *source, 
				const Texture *target,
				NearestNeighborField<Patch, Scalar> *prev,
				typename NearestNeighborField<Patch, Scalar>::Extension *ext, 
				NNSettings &settings) const {
			// assert(source->channels() == target->channels());
			if (source->channels() == channels) {
				return nnf_n<Patch, Scalar, channels>(source, target, prev, ext, settings);
			} else {
				MultiChannelNNF<Patch, Scalar, channels + 1 > op;
				return op(source, target, prev, ext, settings);
			}
		}
	};

#ifndef MAX_SUPPORTED_CHANNELS
#define MAX_SUPPORTED_CHANNELS 12
#endif

	template <typename Patch, typename Scalar>
	struct MultiChannelNNF<Patch, Scalar, MAX_SUPPORTED_CHANNELS + 1 > {

		inline NearestNeighborField<Patch, Scalar> *operator()(
				const Texture *, const Texture *,
				NearestNeighborField<Patch, Scalar> *,
				typename NearestNeighborField<Patch, Scalar>::Extension *,
				NNSettings &) const {
			std::cerr << "Too many channels, we do not support more than " << MAX_SUPPORTED_CHANNELS << " channels!\n";
			return NULL;
		}
	};

	template <typename Patch, typename Scalar>
	NearestNeighborField<Patch, Scalar> *nnf(const Texture *source, const Texture *target,
			NearestNeighborField<Patch, Scalar> *prev,
			typename NearestNeighborField<Patch, Scalar>::Extension *ext, 
			NNSettings &settings) {
		MultiChannelNNF<Patch, Scalar, 1> op;
		return op(source, target, prev, ext, settings);
    }

}

#endif