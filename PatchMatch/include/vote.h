/*******************************************************************************
 * vote.h - voting implementations
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef VOTE_H
#define	VOTE_H

#include <filter.h>
#include <gb.h>
#include <nnf.h>
#include <patch.h>
#include <rng.h>
#include <voting/histogram.h>
#include <voting/meanshift.h>

namespace pm {

#define NOT_NULL(type) reinterpret_cast<type *>(1)

	/**
	 * \brief The method of vote
	 */
	enum VoteMethod {
		/**
		 * \brief Bidirectional Similarity vothing
		 * \see http://www.wisdom.weizmann.ac.il/~vision/VisualSummary.html
		 */
		BiDirSimVoting,
		/**
		 * \brief Bidirectional Similarity with Histogram voting
		 */
		BiDirSimWithHistogramVoting,
		/**
		 * \brief Default voting per pixel (over the corresponding overlapping patches)
		 */
		DefaultVoting,
		/**
		 * \brief Default voting per pixel, with patch gain+bias adjustment
		 */
		DefaultGBAVoting,
		/**
		 * \brief Direct pixel voting
		 */
		DirectVoting,
		/**
		 * \brief Feature-specific voting with histograms for the rest
		 */
		FeatureWithHistogramVoting,
		/**
		 * \brief Voting per pixel, with histogram reweighting
		 */
		HistogramVoting,
		/**
		 * \brief Mean-shift voting per pixel (over the corresponding overlapping patches)
		 */
		MeanShiftVoting,
		/**
		 * \brief Median voting per pixel (over the corresponding overlapping patches)
		 */
		MedianVoting,
		/**
		 * \brief Tiled voting per patch (over the corresponding overlapped pixels)
		 */
		TiledVoting,
		/**
		 * \brief Default voting using a mask as global weight mask
		 */
		WeightedVoting
	};
	
	/**
	 * \brief Type of output weight
	 */
	enum WeightType {
		NoWeight,
		PixelWeight,
		PixelVariance
	};
	
	typedef const float& (*BinFloatOp)(const float&, const float&);

	/**
	 * Parameters for voting
	 */
	struct VoteParams {
		// General:
		int patchSize;
		VoteMethod method;
		Filter filter;
		// Output data:
		WeightType weightType;
		float *weights;
		// BiDir Simimilarity:
		float bidirSimWeight;
		// Mean-shift:
		voting::MeanShiftParams meanShift;
		// Histogram:
		std::vector<int> histChannels;
		std::vector<int> histBins;
		std::vector<Range> histRanges;
		std::vector<float> histWeights;
        std::vector<float> histBoosts;
        bool histNormalize;
		// Weighted:
		Mask weightMask;
        float weightBase;
		// Feature-specific:
		BinFloatOp featureOp;
		float featureDefault;
        std::vector<int> binaryChannels;

		VoteParams() : 
        patchSize(7), 
        method(DefaultVoting), 
        filter(1), 
		weightType(NoWeight),
        weights(NULL),
		bidirSimWeight(0.5f),
		histChannels(),
		histBins(),
		histRanges(),
		histWeights(),
        histBoosts(),
        histNormalize(true),
        weightBase(2.0f),
		featureOp(NULL),
        binaryChannels(){
		}
	};

	template <typename Patch, typename Scalar>
	Image vote(const NearestNeighborField<Patch, Scalar> *nnf, VoteParams &params);

	template <int channels, typename Patch, typename Scalar>
	inline float pixelvar(const NearestNeighborField<Patch, Scalar> *nnf,
		int y, int startY, int endY, 
		int x, int startX, int endX,
		Vec<Scalar, channels> mean) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *target = nnf->target;
		
		// variance sum
		float w = 0;
		int N = 0;
		// for each pixel from the corresponding patches
		for (int py = startY; py < endY; ++py) {
			for (int px = startX; px < endX; ++px) {
				const Patch &patch = nnf->get(py, px);
				int by = y - py; //< pixelInPatch = pos - patchPos
				int bx = x - px;
				PixVal value = target->at<PixVal>(patch.transform(by, bx));
				PixVal diff = mean - value;
				w += diff.dot(diff);
				++N;
			}
		}
		return w / N;
	}
	
	/**
	 * \brief Vote by using a simple filter over overlapping patches
	 * 
	 * \param nnf
	 *			the nearest neighbor field to vote with
	 * \param params
	 *			the voting parameters
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_filter_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *source = nnf->source;
		const Texture *target = nnf->target;
		Image vote = Image::zeros(source->rows, source->cols, IM_32FC(channels));

		const Filter &filter = params.filter;
		// for each pixel of the source
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < source->rows; ++y) {
			for (int x = 0; x < source->cols; ++x) {
                PixVal &votedPixel = vote.at<PixVal>(y, x);
				float w = 0.0f;
				// for each pixel from the corresponding patches
				int startX = std::max(0, x - Patch::width() + 1);
				int endX = std::min(x + 1, nnf->width);
				int startY = std::max(0, y - Patch::width() + 1);
				int endY = std::min(y + 1, nnf->height);
				for (int py = startY; py < endY; ++py) {
					for (int px = startX; px < endX; ++px) {
						const Patch &patch = nnf->get(py, px);
						// /!\ Convolution: the filter is reversed in time
						// => at pos (x,y), the filter is at (0, 0)
						// => at pos (x-dx, y-dy), the filter is at (dx, dy)
						// ++ we assume that only patches on the left / top have a contribution
						int by = y - py; //< pixelInPatch = pos - patchPos
						int bx = x - px;
						PixVal value = target->at<PixVal>(patch.transform(by, bx));
						votedPixel += value * filter[by][bx];
						w += filter[by][bx];
					}
				}
				if (w > 1e-8) {
					votedPixel *= 1.0 / w;
				} else {
#if _OPENMP
#pragma omp critical
#endif
					{
						std::cerr << "w0 @" << y << "/" << x << ", from: " << startY << "/" << startX;
						std::cerr << " to " << endY << "/" << endX << "\n";
					}
				}

				if (votedPixel[0] > 1e10) {
#if _OPENMP
#pragma omp critical
#endif
					{
						std::cerr << "Corrupted vote p@" << y << "/" << x << ": ";
						for (int i = 0; i < channels; ++i) std::cerr << votedPixel[i] << ", ";
						std::cerr << "\n";
					}
				}
                
                // binary channels
                for(int i = 0; i < params.binaryChannels.size(); ++i) {
                    Scalar &v = votedPixel[params.binaryChannels[i]];
                    if(v >= 50.0){
                        v = 100.0;
                    } else {
                        v = 0.0;
                    }
                }

				// store in output weight map if there is one
				switch(params.weightType) {
					case PixelWeight:
						params.weights[source->cols * y + x] = w;
						break;
					case PixelVariance:
						params.weights[source->cols * y + x] = pixelvar(
								nnf,
								y, startY, endY, 
								x, startX, endX, 
								votedPixel
						);
						break;
				}
			}
		}
		return vote;
	}

	/**
	 * \brief Vote by using a simple filter over overlapping patches
	 * and gain+bias adjustment of patches
	 * 
	 * \param nnf
	 *			the nearest neighbor field to vote with
	 * \param params
	 *			the voting parameters
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_filter_gba_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *source = nnf->source;
		const Texture *target = nnf->target;
		Image vote = Image::zeros(source->rows, source->cols, IM_32FC(channels));

		int total = source->rows * source->cols;
		const Filter &filter = params.filter;

		// 1 = compute gain+bias of each target patch
		typedef GainBias<Scalar> GB;
		typedef typename GB::Vec3 Vec3;
		Vec3 *gains = new Vec3[nnf->size], *biases = new Vec3[nnf->size];
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < nnf->height; ++y) {
			for (int x = 0; x < nnf->width; ++x) {
				const typename Patch::OriginalPatchType srcPatch(y, x);
				const Patch &patch = nnf->get(y, x);
				GB::template compute<channels>(source, target, srcPatch, patch,
						gains[nnf->width * y + x], biases[nnf->width * y + x]);
			}
		}

		// for each pixel of the source
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < source->rows; ++y) {
			for (int x = 0; x < source->cols; ++x) {
				float w = 0.0f;
				// for each pixel from the corresponding patches
				int startX = std::max(0, x - Patch::width() + 1);
				int endX = std::min(x + 1, nnf->width);
				int startY = std::max(0, y - Patch::width() + 1);
				int endY = std::min(y + 1, nnf->height);
				for (int py = startY; py < endY; ++py) {
					for (int px = startX; px < endX; ++px) {
						const Patch &patch = nnf->get(py, px);
						// /!\ Convolution: the filter is reversed in time
						// => at pos (x,y), the filter is at (0, 0)
						// => at pos (x-dx, y-dy), the filter is at (dx, dy)
						// ++ we assume that only patches on the left / top have a contribution
						int by = y - py; //< pixelInPatch = pos - patchPos
						int bx = x - px;
						PixVal value = target->at<PixVal>(patch.transform(by, bx));
						// apply gain+bias adjustment
						GB::template applyOn<channels>(value, gains[nnf->width * py + px], biases[nnf->width * py + px]);
						// register the pixel vote
						vote.at<PixVal>(y, x) += value * filter[by][bx];
						w += filter[by][bx];
					}
				}
				if (w > 1e-8) {
					vote.at<PixVal>(y, x) *= 1.0 / w;
				}
				// store in output weight map if there is one
				if (params.weightType == PixelWeight) params.weights[source->cols * y + x] = w;
			}
		}
		return vote;
	}

/**
	 * \brief Direct vote
	 *
	 * \param nnf
	 *			the nearest neighbor field to vote with
	 * \param params
	 *			the voting parameters
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_direct_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *source = nnf->source;
		const Texture *target = nnf->target;
		Image vote = Image::zeros(source->rows, source->cols, IM_32FC(channels));

		int total = source->rows * source->cols;
		if (params.weightType != NoWeight) {
			// store in output weight map if there is one
			switch(params.weightType) {
				case PixelWeight:
					std::fill(params.weights, params.weights + total, 1.0f);
					break;
				case PixelVariance:
					std::fill(params.weights, params.weights + total, 0.0f);
					break;
			}
			
		}
		// for each pixel of the source
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < source->rows; ++y) {
			for (int x = 0; x < source->cols; ++x) {
				// patch position
				int py = std::min(y, nnf->height - 1);
				int px = std::min(x, nnf->width - 1);
				const Patch &patch = nnf->get(py, px);
				// pixel within patch
				int by = y - py;
				int bx = x - px;
				vote.at<PixVal>(y, x) = target->at<PixVal>(patch.transform(by, bx));
			}
		}
		return vote;
	}

	/**
	 * \brief Vote by using a simple filter over overlapping patches
	 * and a tiled strategy with openmp
	 * 
	 * \param nnf
	 *			the nearest neighbor field to vote with
	 * \param params
	 *			the voting parameters
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_filter_tiled_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *source = nnf->source;
		const Texture *target = nnf->target;
		const int h = nnf->height, w = nnf->width; //< nnf bounds
		const int sh = source->rows, sw = source->cols; //< voted image bounds
		Image vote = Image::zeros(sh, sw, IM_32FC(channels));
		int total = source->rows * source->cols;
		float *weights = new float[total](); // /!\ Use default constructor!
		if (params.weights != NULL) {
			switch(params.weightType) {
				case PixelWeight:
					params.weights = weights;
					break;
				case PixelVariance:
					std::fill(params.weights, params.weights + total, 0.0f); //< no valid
					break;
			}
		}
		const Filter &filter = params.filter;
#if _OPENMP
#pragma omp parallel
#endif
		{
			int tid = omp_get_thread_num(), nt = omp_get_num_threads();
			int patchY = std::ceil(float(h) / Patch::width());
			int patchX = std::ceil(float(w) / Patch::width());
			int tiles = patchY * patchX;
			int minTile = tiles * tid / nt;
			int maxTile = tiles * (tid + 1) / nt;
			// for each tile
			for (int tile = minTile; tile < maxTile; ++tile) {
				// for each patch (pixel) of the corresponding tile
				int x = (tile % patchX) * Patch::width();
				int nX = std::min(x + Patch::width(), w);
				for (; x < nX; ++x) {
					int y = (tile / patchX) * Patch::width();
					int nY = std::min(y + Patch::width(), h);
					for (; y < nY; ++y) {
						// /!\ This ^ was done over the nnf field that is smaller than the source!
						const Patch &patch = nnf->get(y, x);
						// for each pixel in the corresponding patch
						for (int px = 0; px < Patch::width(); ++px) {
							for (int py = 0; py < Patch::width(); ++py) {
								// /!\ <-- this is done over the source (nff + patch size)
								PixVal value = target->at<PixVal>(patch.transform(py, px));
								vote.at<PixVal>(y + py, x + px) += value * filter[py][px];
								weights[sw * (y + py) + x + px] += filter[py][px];
							}
						}
					}
				}
			}
		}

		std::cout << ". normalizing\n";
		// normalization of the pixel values
#if _OPENMP
#pragma omp for collapse(2)
#endif
		for (int y = 0; y < sh; ++y) {
			for (int x = 0; x < sw; ++x) {
				float &w = weights[sw * y + x];
				if (w > 1e-8) {
					vote.at<PixVal>(y, x) *= 1.0 / w;
				} else {
					w = 0.0f;
					vote.at<PixVal>(y, x) *= 0.0;
				}
			}
		}
		// free and return
		if (params.weights == NULL) {
			// nobody will still use it => free it!
			delete[] weights;
		}
		return vote;
	}

	template <int channels, typename Patch, typename Scalar>
	Image vote_meanshift_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		typedef voting::Cluster<Scalar, channels> Cluster;
		const Texture *source = nnf->source;
		const Texture *target = nnf->target;
		Image vote = Image::zeros(source->rows, source->cols, IM_32FC(channels));
		// for each pixel of the source
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < source->rows; ++y) {
			for (int x = 0; x < source->cols; ++x) {
				// for each pixel from the corresponding patches
				int startX = std::max(0, x - Patch::width() + 1);
				int endX = std::min(x + 1, nnf->width);
				int startY = std::max(0, y - Patch::width() + 1);
				int endY = std::min(y + 1, nnf->height);
				// mean-shift voting on the feature space
				// generated by the pixels on (startY:endY-1, startX:endX-1)
				// - use params.meanShiftWindows to discard
				typedef std::vector<PixVal> PixCluster;
				// workspace
				int pw = endX - startX;
				int ph = endY - startY;
				int N = pw * ph;
				std::vector<PixVal> points(N, PixVal()); //< data points
				int pidx = 0;
				for (int py = startY; py < endY; ++py) {
					for (int px = startX; px < endX; ++px, ++pidx) {
						const Patch &patch = nnf->get(py, px);
						points[pidx] = target->at<PixVal>(patch.transform(y - py, x - px)); // * filter[by][bx];
					}
				}
				// mean-shift algorithm
				std::vector<Cluster> clusters;
				clusters.reserve(5);
				voting::meanshift(points, clusters, params.meanShift);

				// pixel vote using the best cluster
				int best = 0;
				for (int i = 1, n = clusters.size(); i < n; ++i) {
					if (clusters[i].votes > clusters[best].votes) {
						best = i;
					}
				}
				vote.at<PixVal>(y, x) = clusters[best].mean;

				switch(params.weightType) {
					case PixelWeight:
						params.weights[vote.cols * y + x] = float(clusters[best].votes);
						break;
					case PixelVariance:
						params.weights[vote.cols * y + x] = pixelvar(
								nnf,
								y, startY, endY, 
								x, startX, endX, 
								vote.at<PixVal>(y, x)
						);
						break;
				}
			}
		}
		return vote;
	}
	
	template <int K, typename Scalar>
	class MedianList{
	public:
		typedef float Weight;
		typedef std::pair<Scalar, Weight> Item;
		struct ItemOrder {
			bool operator()(const Item &left, const Item &right) {
				return left.first < right.first;
			}
		};
		typedef Vec<Scalar, K> Value;
		typedef Vec<float, K> WeightVec;
		MedianList(int N) : count(0), total(0.0f) {
			for(int k = 0; k < K; ++k) lists[k].resize(N); // allocate space
		}
		inline void push(const Value &v, Weight w) {
			for(int k = 0; k < K; ++k){
				lists[k][count] = std::make_pair(v[k], w);
			}
			++count;
			total += w;
		}
		inline void sort() {
			for(int k = 0; k < K; ++k) std::sort(lists[k].begin(), lists[k].end(), ItemOrder());
		}
		inline Value get() const {
			const Weight middle = total * 0.5f; // varies with the position (since boundary patches have smaller weights)
			Value val; // value computation
			WeightVec weight = WeightVec::zeros(); // weight storage
			for(int i = 0; i < count; ++i) {
				// separately for each channel
				bool more = false;
				for(int k = 0; k < K; ++k) {
					const Item &it = lists[k][i];
					if(weight[k] < middle){
						weight[k] += it.second;
						// if we went past the center, that's the median
						if(weight[k] >= middle) val[k] = it.first;
						else more = true; // else we need to go farther
					} // else nothing more
				}
				if(!more) break;
			}
			return val;
		}
	private:
		int count;
		Weight total;
		std::vector<Item> lists[K];
	};
	
	template <int channels, typename Patch, typename Scalar>
	Image vote_median_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *source = nnf->source;
		const Texture *target = nnf->target;
		Image vote = Image::zeros(source->rows, source->cols, IM_32FC(channels));

		// weights recording
		if (params.weights != NULL) {
			std::fill(params.weights, params.weights + (source->cols * source->rows), 0.0f);
		}
		
		// for each pixel of the source
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < source->rows; ++y) {
			for (int x = 0; x < source->cols; ++x) {
				// for each pixel from the corresponding patches
				int startX = std::max(0, x - Patch::width() + 1);
				int endX = std::min(x + 1, nnf->width);
				int startY = std::max(0, y - Patch::width() + 1);
				int endY = std::min(y + 1, nnf->height);
				// sort pixel values
				int pw = endX - startX;
				int ph = endY - startY;
				int N = pw * ph;
				MedianList<channels, Scalar> points(N); //< data points
				int pidx = 0;
				for (int py = startY; py < endY; ++py) {
					for (int px = startX; px < endX; ++px, ++pidx) {
						const Patch &patch = nnf->get(py, px);
						int by = y - py; // location in the patch
						int bx = x - px;
						points.push(target->at<PixVal>(patch.transform(by, bx)), params.filter[by][bx]);
					}
				}
				// sort data, channel by channel
				points.sort();
				// get the sort-of median
				vote.at<PixVal>(y, x) = points.get();
				
				switch(params.weightType) {
					case PixelWeight:
						params.weights[vote.cols * y + x] = 1.0f;
						break;
					case PixelVariance:
						params.weights[vote.cols * y + x] = pixelvar(
								nnf,
								y, startY, endY, 
								x, startX, endX, 
								vote.at<PixVal>(y, x)
						);
						break;
				}
			}
		}
		return vote;
	}

	/**
	 * \brief Vote by using a weighted filter over overlapping patches
	 * that is reweighted as the histogram distribution changes.
	 * 
	 * \param nnf
	 *			the nearest neighbor field to vote with
	 * \param params
	 *			the voting parameters
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_histogram_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *T = nnf->source;
		const Texture *E = nnf->target;
		Image vote = Image::zeros(T->rows, T->cols, IM_32FC(channels));

		int numPixels = T->rows * T->cols;
        float areaT = numPixels, areaE = E->rows * E->cols;
		const Filter &filter = params.filter;

		// 1 = compute the histograms
		typedef voting::Histogram<Patch, Scalar> Hist;
		const int K = params.histChannels.size();
		std::vector<Hist> histT(K); // to be changed
		std::vector<Hist> histE(K); // DO NOT change
		for (int i = 0; i < K; ++i) {
			histT[i] = Hist(T, params.histChannels[i], params.histBins[i], params.histRanges[i]);
			histE[i] = Hist(E, params.histChannels[i], params.histBins[i], params.histRanges[i]);
		}

		// 2 = random traversal
		std::vector<Point<int> > index;
		index.resize(numPixels);
		for (int y = 0, i = 0; y < T->rows; ++y) {
			for (int x = 0; x < T->cols; ++x, ++i) {
				index[i] = Point2i(x, y);
			}
		}
		knuth_shuffle(unif01, &index[0], numPixels);

		// 3 = histogram-weighted voting
		// for each pixel of the source
		float ratioTE = areaT / areaE;
		for (int idx = 0; idx < numPixels; ++idx) {
			// the index
			int y = index[idx].y;
			int x = index[idx].x;
			// std::cout << "#" << idx << " @(" << y << " / " << x << ")\n";
			float weightSum = 0.0f;
			PixVal &votedPixel = vote.at<PixVal>(y, x);
			// for each pixel from the corresponding patches
			int startX = std::max(0, x - Patch::width() + 1);
			int endX = std::min(x + 1, nnf->width);
			int startY = std::max(0, y - Patch::width() + 1);
			int endY = std::min(y + 1, nnf->height);
			for (int py = startY; py < endY; ++py) {
				for (int px = startX; px < endX; ++px) {
					const Patch &patch = nnf->get(py, px);
					int by = y - py; //< pixelInPatch = pos - patchPos
					int bx = x - px;
					PixVal value = E->at<PixVal>(patch.transform(by, bx));
					// reweighting using the histogram
                    // w(p) = w * (1 + \sum_k a_k [hb_trg - hb_src]+) / (1 + \sum_k b_k [hb_src - hb_trg]+)
					float dA = 1.0f, dB = 1.0f;
					for (int k = 0; k < K; ++k) {
						int ch = params.histChannels[k];
                        int hbT = histT[k].count(value[ch]);
                        int hbE = histE[k].count(value[ch]);
                        float da;
                        if (params.histNormalize) {
                            da = float(hbE) * ratioTE - float(hbT);
                        } else {
                            da = hbE - hbT;
                        }
                        float db = -da;
                        if (params.histBoosts.size() > 0 && da > 0) dA += params.histBoosts[k] * da;
						if (db > 0) dB += params.histWeights[k] * db;
					}
					float w = filter[by][bx] * dA / dB;
					votedPixel += value * w;
					weightSum += w;
				}
			}
			if (weightSum > 1e-8) {
				votedPixel *= 1.0 / weightSum;
			} else {
				std::cerr << "Dangerously low weight @" << y << "/" << x << " = " << weightSum << "\n";
			}
            
            // binary channels
            for(int i = 0; i < params.binaryChannels.size(); ++i) {
                Scalar &v = votedPixel[params.binaryChannels[i]];
                if(v >= 50.0){
                    v = 100.0;
                } else {
                    v = 0.0;
                }
            }

			// update the source histograms
			const Scalar *prevPixel = T->ptr<Scalar>(y, x);
			for (int k = 0; k < K; ++k) {
				Hist &h = histT[k];
				int ch = params.histChannels[k];
				// previous value
				float v = prevPixel[ch];
				h.count(v) -= 1;
				// new value
				v = votedPixel[ch];
				h.count(v) += 1;
			}

			// store in output weight map if there is one
			switch(params.weightType) {
                case PixelWeight:
                    params.weights[vote.cols * y + x] = weightSum;
                    break;
                case PixelVariance:
                    params.weights[vote.cols * y + x] = pixelvar(
						nnf,
						y, startY, endY, 
						x, startX, endX, 
						votedPixel
                    );
                    break;
			}
		}
		return vote;
	}
	
	/**
	 * \brief Vote using a mask as weight for the pixels
	 * 
	 * \param nnf
	 *			the nearest neighbor field to vote with
	 * \param params
	 *			the voting parameters
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_weighted_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *source = nnf->source;
		const Texture *target = nnf->target;
		Image vote = Image::zeros(source->rows, source->cols, IM_32FC(channels));

		const Mask &mask = params.weightMask;
        const float base = params.weightBase;
        const int midP = Patch::width() / 2;
		// for each pixel of the source
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < source->rows; ++y) {
			for (int x = 0; x < source->cols; ++x) {
                PixVal &votedPixel = vote.at<PixVal>(y, x);
				float wSum = 0.0f;
				// pixel distance
				float pxDist = mask.at<float>(
					std::min(std::max(y - midP, 0), mask.rows - 1), 
					std::min(std::max(x - midP, 0), mask.cols - 1)
				);
				// for each pixel from the corresponding patches
				int startX = std::max(0, x - Patch::width() + 1);
				int endX = std::min(x + 1, nnf->width);
				int startY = std::max(0, y - Patch::width() + 1);
				int endY = std::min(y + 1, nnf->height);
				for (int py = startY; py < endY; ++py) {
					for (int px = startX; px < endX; ++px) {
						const Patch &patch = nnf->get(py, px);
						// /!\ Convolution: the filter is reversed in time
						// => at pos (x,y), the filter is at (0, 0)
						// => at pos (x-dx, y-dy), the filter is at (dx, dy)
						// ++ we assume that only patches on the left / top have a contribution
						int by = y - py; //< pixelInPatch = pos - patchPos
						int bx = x - px;
						PixVal value = target->at<PixVal>(patch.transform(by, bx));
						float w = std::pow(base, -(mask.at<float>(py, px) - pxDist));
						votedPixel += value * w;
						wSum += w;
					}
				}
				if (wSum > 1e-8) {
					votedPixel *= 1.0 / wSum;
				} else {
#if _OPENMP
#pragma omp critical
#endif
					{
						std::cerr << "w0 (" << wSum << ") @" << y << "/" << x << ", from: " << startY << "/" << startX;
						std::cerr << " to " << endY << "/" << endX << "\n";
                        for (int py = startY; py < endY; ++py) {
                            for (int px = startX; px < endX; ++px) {
                                float w = std::pow(2.0f, -(mask.at<float>(py, px) - pxDist));
                                std::cerr << "(" << py << ", " << px << ") -> " << mask.at<float>(py, px) << " -> " << w << ", ";
                            }
                        }
                        std::cerr << "\n\n";
					}
				}
                
                // binary channels
                for(int i = 0; i < params.binaryChannels.size(); ++i) {
                    Scalar &v = votedPixel[params.binaryChannels[i]];
                    if(v >= 50.0){
                        v = 100.0;
                    } else {
                        v = 0.0;
                    }
                }

				// store in output weight map if there is one
				switch(params.weightType) {
					case PixelWeight:
						params.weights[vote.cols * y + x] = wSum;
						break;
					case PixelVariance:
						params.weights[vote.cols * y + x] = pixelvar(
								nnf,
								y, startY, endY, 
								x, startX, endX, 
								votedPixel
						);
						break;
				}
			}
		}
		return vote;
	}
	
	/**
	 * \brief Vote using a bidirectional similarity measure
	 * 
	 * \param nnfs
	 *			array containing [0]:the T to S nnf, [1]: the S to T nnf
	 * \param params
	 *			the voting parameters including the S to T nnf
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_bidir_sim_n(const NearestNeighborField<Patch, Scalar> *nnfs,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		typedef NearestNeighborField<Patch, Scalar> NNF;
		typedef typename NNF::SourcePatch SourcePatch;
		typedef typename NNF::TargetPatch TargetPatch;
		const NNF &dirNNF = nnfs[0]; // direct T to S nnf
		const NNF &revNNF = nnfs[1]; // reverse S to T nnf
		const Texture *T = dirNNF.source; // T in bidir sim
		const Texture *S = dirNNF.target; // S in bidir sim
		Image vote = Image::zeros(T->rows, T->cols, IM_32FC(channels));
		const Filter &filter = params.filter;
		float *weights;
		// set weights to 0 to start
		if (params.weights != NULL) {
			weights = params.weights;
			std::fill(weights, weights + (T->height * T->width), 0.0f);
		} else {
			weights = new float[T->height * T->width](); // init to 0
		}
		
		// group weights
		float dirWeight, revWeight;
		{
			double regularizer = std::max(dirNNF.size, revNNF.size);
			double dw = (1.0 - params.bidirSimWeight) * regularizer / dirNNF.size;
			double rw = params.bidirSimWeight * regularizer / revNNF.size;
			dirWeight = float(dw);
			revWeight = float(rw);
			std::cout << "dw: " << dirWeight << ", rw: " << revWeight << "\n";
		}
		
		// process the reverse nnf first
		for (int y = 0; y < revNNF.height; ++y) {
			for(int x = 0; x < revNNF.width; ++x) {
				SourcePatch p(y, x);
				const TargetPatch &q = revNNF.get(y, x);
				// for each pixel from the patch p in S
				for (typename SourcePatch::IndexIterator it = p.begin(); it; ++it) {
					typename SourcePatch::Index i = *it;
					const Point2i lp(p.transform(i));
					PixVal pix = S->at<PixVal>(lp); // transformation in S space
					// target in T
					const Point2i lq(q * i); // q.transform(i)
					const Point2i delta = lp - Point2i(x, y); // (x,y) is at the top-left of the patch
					float w = filter[delta.y][delta.x] * revWeight;
					vote.at<PixVal>(lq.y, lq.x) += pix * w;
					// store weight (count for Ns when filter = 1)
					weights[vote.cols * lq.y + lq.x] += w;
				}
			}
		}
		
		// process the direct nnf finally
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < T->rows; ++y) {
			for (int x = 0; x < T->cols; ++x) {
                PixVal &votedPixel = vote.at<PixVal>(y, x);
				float wT = 0.0f;
				// for each pixel from the corresponding patches
				int startX = std::max(0, x - Patch::width() + 1);
				int endX = std::min(x + 1, dirNNF.width);
				int startY = std::max(0, y - Patch::width() + 1);
				int endY = std::min(y + 1, dirNNF.height);
				for (int py = startY; py < endY; ++py) {
					for (int px = startX; px < endX; ++px) {
						const TargetPatch &patch = dirNNF.get(py, px);
						// /!\ Convolution: the filter is reversed in time
						// => at pos (x,y), the filter is at (0, 0)
						// => at pos (x-dx, y-dy), the filter is at (dx, dy)
						// ++ we assume that only patches on the left / top have a contribution
						int by = y - py; //< pixelInPatch = pos - patchPos
						int bx = x - px;
						PixVal value = S->at<PixVal>(patch.transform(by, bx));
						float w = filter[by][bx] * dirWeight;
						votedPixel += value * w;
						wT += w;
					}
				}
				float wS = weights[T->cols * y + x];
				float wSum = wT + wS;
				weights[T->cols * y + x] = wSum;
				if (wSum > 1e-8) {
					votedPixel *= 1.0 / wSum;
				} else {
#if _OPENMP
#pragma omp critical
#endif
					{
						std::cerr << "w0 @" << y << "/" << x << ", from: " << startY << "/" << startX;
						std::cerr << " to " << endY << "/" << endX << "\n";
					}
				}
                
                // binary channels
                for(int i = 0; i < params.binaryChannels.size(); ++i) {
                    Scalar &v = votedPixel[params.binaryChannels[i]];
                    if(v >= 50.0){
                        v = 100.0;
                    } else {
                        v = 0.0;
                    }
                }
			}
		}
		
		// free weights if only used locally
		if(params.weights == NULL) {
			delete[] weights;
		}
		return vote;
	}
	
	/**
	 * \brief Vote using a bidirectional similarity measure and histograms
	 * 
	 * \param nnfs
	 *			array containing [0]:the T to S nnf, [1]: the S to T nnf
	 * \param params
	 *			the voting parameters including the S to T nnf
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_bidir_sim_histogram_n(const NearestNeighborField<Patch, Scalar> *nnfs,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		typedef NearestNeighborField<Patch, Scalar> NNF;
		typedef typename NNF::SourcePatch SourcePatch;
		typedef typename NNF::TargetPatch TargetPatch;
		const NNF &dirNNF = nnfs[0]; // direct T to S nnf
		const NNF &revNNF = nnfs[1]; // reverse S to T nnf
		const Texture *T = dirNNF.source; // T in bidir sim
		const Texture *E = dirNNF.target; // S in bidir sim
		Image vote = Image::zeros(T->rows, T->cols, IM_32FC(channels));
		// const Filter &filter = params.filter;
		
		float *weights;
		// set weights to 0 to start
		if (params.weights != NULL) {
			weights = params.weights;
			std::fill(weights, weights + (T->height * T->width), 0.0f);
		} else {
			weights = new float[T->height * T->width](); // init to 0
		}
		
		// 1 = constructor accumulator from revNNF
		std::vector<Point2i> *accum = new std::vector<Point2i>[T->height * T->width](); // init vectors
		if(accum == NULL){
			mexErrMsgIdAndTxt("MATLAB:vote:new", "Failed allocation of accum in BidirSimVoting with Histograms");
		}
		// using the reverse nnf first
		for (int y = 0; y < revNNF.height; ++y) {
			for(int x = 0; x < revNNF.width; ++x) {
				SourcePatch p(y, x);
				const TargetPatch &q = revNNF.get(y, x);
				// for each pixel from the patch p in S
				for (typename SourcePatch::IndexIterator it = p.begin(); it; ++it) {
					typename SourcePatch::Index i = *it;
					// pixel location in E
					const Point2i lp(p.transform(i));
					// pixel location in T
					const Point2i lq(q.transform(i));
					accum[vote.cols * lq.y + lq.x].push_back(lp);
				}
			}
		}
		
		// group weights
		float dirWeight, revWeight;
		{
			double regularizer = std::max(dirNNF.size, revNNF.size);
			double dw = (1.0 - params.bidirSimWeight) * regularizer / dirNNF.size;
			double rw = params.bidirSimWeight * regularizer / revNNF.size;
			dirWeight = float(dw);
			revWeight = float(rw);
			std::cout << "dw: " << dirWeight << ", rw: " << revWeight << "\n";
		}
		
		// histogram voting now, using both the direct NNF and the accumulator
		int numPixels = T->rows * T->cols;
        float areaT = numPixels, areaE = E->rows * E->cols;
		float ratioTE = areaT / areaE;
		
		// 2 = compute the histograms
		typedef voting::Histogram<Patch, Scalar> Hist;
		const int K = params.histChannels.size();
		std::vector<Hist> histT(K); // to be changed
		std::vector<Hist> histE(K); // DO NOT change
		for (int i = 0; i < K; ++i) {
			histT[i] = Hist(T, params.histChannels[i], params.histBins[i], params.histRanges[i]);
			histE[i] = Hist(E, params.histChannels[i], params.histBins[i], params.histRanges[i]);
		}

		// 2 = random traversal
		std::vector<Point2i> index;
		index.resize(numPixels);
		for (int y = 0, i = 0; y < T->rows; ++y) {
			for (int x = 0; x < T->cols; ++x, ++i) {
				index[i] = Point2i(x, y);
			}
		}
		knuth_shuffle(unif01, &index[0], numPixels);

		// 3 = histogram-weighted voting
		// for each pixel of the source
		//  + the accumulated pixels with revNNF
		for (int idx = 0; idx < numPixels; ++idx) {
			// the index
			int y = index[idx].y;
			int x = index[idx].x;
			float weightSum = 0.0f;
			PixVal &votedPixel = vote.at<PixVal>(y, x);
			// for each pixel from the corresponding patches
			int startX = std::max(0, x - Patch::width() + 1);
			int endX = std::min(x + 1, dirNNF.width);
			int startY = std::max(0, y - Patch::width() + 1);
			int endY = std::min(y + 1, dirNNF.height);
			for (int py = startY; py < endY; ++py) {
				for (int px = startX; px < endX; ++px) {
					const Patch &patch = dirNNF.get(py, px);
					int by = y - py; //< pixelInPatch = pos - patchPos
					int bx = x - px;
					PixVal value = E->at<PixVal>(patch.transform(by, bx));
					
					// reweighting using the histogram
                    // w(p) = w * (1 + \sum_k a_k [hb_trg - hb_src]+) / (1 + \sum_k b_k [hb_src - hb_trg]+)
					float dH = 1.0f;
					for (int k = 0; k < K; ++k) {
						int ch = params.histChannels[k];
                        int hbT = histT[k].count(value[ch]);
                        int hbE = histE[k].count(value[ch]);
                        float dh;
                        if (params.histNormalize) {
                            dh = float(hbT) - float(hbE) * ratioTE; // relatively to T
                        } else {
                            dh = hbT - hbE;
                        }
						if (dh > 0) dH += params.histWeights[k] * dh;
					}
					// filter[by][bx] / dH
					float w = dirWeight / dH; // Note: filter is not used in E!
					votedPixel += value * w;
					weightSum += w;
				}
			}
			// for each pixel from the accumulated map
			const std::vector<Point2i> &pixelAccum = accum[vote.cols * y + x];
			for (int i = 0, n = pixelAccum.size(); i < n; ++i) {
				const Point2i &p = pixelAccum[i];
				PixVal value = E->at<PixVal>(p);

				// reweighting using the histogram
				// w(p) = w * (1 + \sum_k a_k [hb_trg - hb_src]+) / (1 + \sum_k b_k [hb_src - hb_trg]+)
				float dH = 1.0f;
				for (int k = 0; k < K; ++k) {
					int ch = params.histChannels[k];
					int hbT = histT[k].count(value[ch]);
					int hbE = histE[k].count(value[ch]);
					float dh;
					if (params.histNormalize) {
						dh = float(hbT) - float(hbE) * ratioTE; // relatively to T
					} else {
						dh = hbE - hbT;
					}
					if (dh > 0) dH += params.histWeights[k] * dh;
				}
				float w = revWeight / dH;
				votedPixel += value * w;
				weightSum += w;
			}
			
			if (weightSum > 1e-8) {
				votedPixel *= 1.0 / weightSum;
			} else {
				std::cerr << "Dangerously low weight @" << y << "/" << x << " = " << weightSum << "\n";
			}
            
            // binary channels
            for(int i = 0; i < params.binaryChannels.size(); ++i) {
                Scalar &v = votedPixel[params.binaryChannels[i]];
                if(v >= 50.0){
                    v = 100.0;
                } else {
                    v = 0.0;
                }
            }

			// update the source histograms
			const Scalar *prevPixel = T->ptr<Scalar>(y, x);
			for (int k = 0; k < K; ++k) {
				float v;
				Hist &h = histT[k];
				int ch = params.histChannels[k];
				// previous value
				v = prevPixel[ch];
				h.count(v) -= 1;
				// new value
				v = votedPixel[ch];
				h.count(v) += 1;
			}

			// store in output weight map if there is one
			switch(params.weightType) {
                case PixelWeight:
                    params.weights[vote.cols * y + x] = weightSum;
                    break;
                case PixelVariance:
                    params.weights[vote.cols * y + x] = pixelvar(
                            &dirNNF,
                            y, startY, endY, 
                            x, startX, endX, 
                            votedPixel
                    );
                    break;
			}
		}
		// free memory
		delete[] accum;
		if(params.weights == NULL) {
			delete[] weights;
		}
		
		// done!
		return vote;
	}
	
	/**
	 * \brief Feature-specific voting and histogram voting for the rest
	 * 
	 * \param nnf
	 *			the nearest neighbor field to vote with
	 * \param params
	 *			the voting parameters
	 * \return the voted picture
	 */
	template <int channels, typename Patch, typename Scalar>
	Image vote_feature_histogram_n(const NearestNeighborField<Patch, Scalar> *nnf,
			VoteParams &params) {
		typedef Vec<Scalar, channels> PixVal;
		const Texture *T = nnf->source;
		const Texture *E = nnf->target;
		Image vote = Image::zeros(T->rows, T->cols, IM_32FC(channels));

		int numPixels = T->rows * T->cols;
        float areaT = numPixels, areaE = E->rows * E->cols;
		const Filter &filter = params.filter;

		// 1 = compute the histograms
		typedef voting::Histogram<Patch, Scalar> Hist;
		const int K = params.histChannels.size();
		std::vector<Hist> histT(K); // to be changed
		std::vector<Hist> histE(K); // DO NOT change
		for (int i = 0; i < K; ++i) {
			histT[i] = Hist(T, params.histChannels[i], params.histBins[i], params.histRanges[i]);
			histE[i] = Hist(E, params.histChannels[i], params.histBins[i], params.histRanges[i]);
		}
		
		// the feature channels
		std::vector<unsigned int> features;
		for(unsigned int c = 0; c < channels; ++c) {
			bool isFeature = true;
			for(unsigned int k = 0; k < K; ++k) {
				if(c == params.histChannels[k]){
					std::cout << "H" << c << ", ";
					isFeature = false;
					break;
				}
			}
			if(isFeature) {
				std::cout << "F" << c << ", ";
				features.push_back(c);
			}
		}
		std::cout << "\n";

		// 2 = random traversal
		std::vector<Point<int> > index;
		index.resize(numPixels);
		for (int y = 0, i = 0; y < T->rows; ++y) {
			for (int x = 0; x < T->cols; ++x, ++i) {
				index[i] = Point2i(x, y);
			}
		}
		knuth_shuffle(unif01, &index[0], numPixels);

		// 3 = histogram-weighted voting
		// for each pixel of the source
		float ratioTE = areaT / areaE;
		for (int idx = 0; idx < numPixels; ++idx) {
			// the index
			int y = index[idx].y;
			int x = index[idx].x;
			// std::cout << "#" << idx << " @(" << y << " / " << x << ")\n";
			float weightSum = 0.0f;
			PixVal &votedPixel = vote.at<PixVal>(y, x);
			for(int f = 0; f < features.size(); ++f) {
				votedPixel[features[f]] = FLT_MAX;
			}
			// for each pixel from the corresponding patches
			int startX = std::max(0, x - Patch::width() + 1);
			int endX = std::min(x + 1, nnf->width);
			int startY = std::max(0, y - Patch::width() + 1);
			int endY = std::min(y + 1, nnf->height);
			for (int py = startY; py < endY; ++py) {
				for (int px = startX; px < endX; ++px) {
					const Patch &patch = nnf->get(py, px);
					int by = y - py; //< pixelInPatch = pos - patchPos
					int bx = x - px;
					PixVal value = E->at<PixVal>(patch.transform(by, bx));
					// extract feature data
					PixVal F;
					for(int f = 0; f < features.size(); ++f){
						int c = features[f];
						F[c] = value[c];
						value[c] = Scalar(0.0);
					}
					// reweighting using the histogram
                    // w(p) = w * (1 + \sum_k a_k [hb_trg - hb_src]+) / (1 + \sum_k b_k [hb_src - hb_trg]+)
					float dA = 1.0f, dB = 1.0f;
					for (int k = 0; k < K; ++k) {
						int ch = params.histChannels[k];
                        int hbT = histT[k].count(value[ch]);
                        int hbE = histE[k].count(value[ch]);
                        float da;
                        if (params.histNormalize) {
                            da = float(hbE) * ratioTE - float(hbT);
                        } else {
                            da = hbE - hbT;
                        }
                        float db = -da;
                        if (params.histBoosts.size() > 0 && da > 0) dA += params.histBoosts[k] * da;
						if (db > 0) dB += params.histWeights[k] * db;
					}
					float w = filter[by][bx] * dA / dB;
					votedPixel += value * w;
					weightSum += w;
					// feature specific values
					for(int f = 0; f < features.size(); ++f) {
						int c = features[f];
						votedPixel[c] = params.featureOp(votedPixel[c], F[c]);
					}
				}
			}
			float invWeight = 1.0f / weightSum;
			for(int i = 0; i < K; ++i) votedPixel[params.histChannels[i]] *= invWeight; // c cannot intersect with k!
			if (weightSum <= 1e-8) {
				std::cerr << "Dangerously low weight @" << y << "/" << x << " = " << weightSum << "\n";
			}
            // binary channels
            for(int i = 0; i < params.binaryChannels.size(); ++i) {
                Scalar &v = votedPixel[params.binaryChannels[i]];
                if(v >= 50.0){
                    v = 100.0;
                } else {
                    v = 0.0;
                }
            }

			// update the source histograms
			const Scalar *prevPixel = T->ptr<Scalar>(y, x);
			for (int k = 0; k < K; ++k) {
				Hist &h = histT[k];
				int ch = params.histChannels[k];
				// previous value
				float v = prevPixel[ch];
				h.count(v) -= 1;
				// new value
				v = votedPixel[ch];
				h.count(v) += 1;
			}

			// store in output weight map if there is one
			switch(params.weightType) {
                case PixelWeight:
                    params.weights[vote.cols * y + x] = weightSum;
                    break;
                case PixelVariance:
                    params.weights[vote.cols * y + x] = pixelvar(
						nnf,
						y, startY, endY, 
						x, startX, endX, 
						votedPixel
                    );
                    break;
			}
		}
		return vote;
	}

	// #########################################################################
	// ##### Multi-channels voting #############################################
	// #########################################################################

	template <typename Patch, typename Scalar, int channels>
	struct MultiChannelVote {

		inline Image operator()(const NearestNeighborField<Patch, Scalar> *nnf, VoteParams &params) const {
			if (nnf->source->channels() == channels) {
				switch (params.method) {
					case BiDirSimVoting:
						return vote_bidir_sim_n<channels>(nnf, params);
					case BiDirSimWithHistogramVoting:
						return vote_bidir_sim_histogram_n<channels>(nnf, params);
					case DefaultVoting:
						return vote_filter_n<channels>(nnf, params);
					case DefaultGBAVoting:
						return vote_filter_gba_n<channels>(nnf, params);
					case DirectVoting:
						return vote_direct_n<channels>(nnf, params);
					case FeatureWithHistogramVoting:
						return vote_feature_histogram_n<channels>(nnf, params);
					case HistogramVoting:
						return vote_histogram_n<channels>(nnf, params);
					case MeanShiftVoting:
						return vote_meanshift_n<channels>(nnf, params);
					case MedianVoting:
						return vote_median_n<channels>(nnf, params);
					case TiledVoting:
						return vote_filter_tiled_n<channels>(nnf, params);
					case WeightedVoting:
						return vote_weighted_n<channels>(nnf, params);
					default:
						std::cerr << "Invalid voting method: " << params.method << " !";
						return Image();
				}
			} else {
				MultiChannelVote<Patch, Scalar, channels + 1 > vop;
				return vop(nnf, params);
			}
		}
	};

#ifndef MAX_SUPPORTED_CHANNELS
#define MAX_SUPPORTED_CHANNELS 12
#endif

	template <typename Patch, typename Scalar>
	struct MultiChannelVote<Patch, Scalar, MAX_SUPPORTED_CHANNELS + 1 > {

		inline Image operator()(const NearestNeighborField<Patch, Scalar> *nnf, VoteParams &params) const {
			std::cerr << "Warning: voting only supported up to " << MAX_SUPPORTED_CHANNELS << " channels!\n";
			return Image();
		}
	};

	template <typename Patch, typename Scalar>
	Image vote(const NearestNeighborField<Patch, Scalar> *nnf, VoteParams &params) {
		MultiChannelVote<Patch, Scalar, 1> vop;
		return vop(nnf, params);
	}
}

#endif	/* VOTE_H */

