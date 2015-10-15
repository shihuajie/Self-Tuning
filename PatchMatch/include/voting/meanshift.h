/*******************************************************************************
 * meanshift.h - mean-shift algorithm for voting
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef MEANSHIFT_H
#define	MEANSHIFT_H

#include <rng.h>
#include <vector>

namespace pm {

	namespace voting {

		/**
		 * \brief Mean-shift parameters
		 */
		struct MeanShiftParams {
			std::vector<float> windows;
			float threshold;
			bool shuffle;
			// default constructor
			MeanShiftParams() : threshold(1.0f), shuffle(true){}
		};

		template <typename Scalar, int N>
		inline bool inWindow(const Vec<Scalar, N> &mean, const Vec<Scalar, N> &point,
				const std::vector<float> &winSizes) {
			bool res = true;
			for (int i = 0; i < N; ++i) {
				if (std::abs(mean[i] - point[i]) > winSizes[i]) return false;
			}
			return res;
		}

		template <typename Scalar, int N>
		inline Scalar windowDist(const Vec<Scalar, N> &q, const Vec<Scalar, N> p,
				const std::vector<float> &winSizes) {
			Scalar dist = 0;
			for (int i = 0; i < N; ++i) {
				dist += std::abs(q[i] - p[i]) / winSizes[i];
			}
			return dist;
		}

		/**
		 * \brief Single cluster description
		 */
		template <typename Scalar, int N>
		struct Cluster {
			typedef Vec<Scalar, N> PixVal;
			/**
			 * \brief The cluster mean
			 */
			PixVal mean;
			/**
			 * \brief The number of votes for that cluster
			 */
			int votes;
		};

		/**
		 * \brief Mean-shift algorithm computing clusters of points within
		 * a feature space by iteratively updating a cluster mean within
		 * a window of the feature space until convergence
		 * 
		 * \param points
		 *			the feature points
		 * \param clusters
		 *			the cluster storage
		 * \param params
		 *			the algorithm parameters
		 */
		template <typename Scalar, int Channels>
		void meanshift(const std::vector<Vec<Scalar, Channels> > &points,
				std::vector<Cluster<Scalar, Channels> > &clusters,
				const MeanShiftParams &params) {
			typedef Vec<Scalar, Channels> PixVal;
			typedef Cluster<Scalar, Channels> Clu;

			// workspace creation
			const int N = points.size();
			std::vector<bool> visited(N, false); //< visit flags
			std::vector<int> perm(N, 0); //< visit permutation
			for (int i = 0; i < N; ++i) perm[i] = i;
			// do we use a permutation or a deterministic index
			if(params.shuffle){
				knuth_shuffle(pm::unif01, &perm[0], N);
			}

			// we visit at most that many possible clusters
			for (int pidx = 0; pidx < N; ++pidx) {
				if (visited[perm[pidx]]) continue; // we already visited that one once
				int i = perm[pidx], it = 0;
				PixVal lastMean = points[i];
				visited[i] = true;
				// until the mean converges
				PixVal diff;
				int votes;
				do {
					// compute the new mean
					PixVal newMean = PixVal::all(0.0);
					votes = 0;
					for (int k = 0; k < N; ++k) {
						// if a point is within window, we take it
						if (inWindow(lastMean, points[k], params.windows)) {
							newMean += points[k];
							visited[k] = true;
							++votes;
						}
					}
					newMean *= 1 / Scalar(votes);
					// mean displacement
					diff = lastMean - newMean;
					// update the mean
					lastMean = newMean;
				} while (++it < 100 || diff.dot(diff) > params.threshold);
				Clu clu = {
					lastMean,
					votes
				};

				// try to merge it with other clusters
				Scalar cdist = std::numeric_limits<Scalar>().max();
				int cidx = -1;
				for (int c = 0, n = clusters.size(); c < n; ++c) {
					Scalar dist = windowDist(lastMean, clusters[i].mean, params.windows);
					if (dist < cdist) {
						cdist = dist;
						cidx = c;
					}
				}
				if (cdist < Channels) {
					// merge
					clusters[cidx].votes += clu.votes; // only?
					clusters[cidx].mean = (clu.mean + clusters[cidx].mean) * 0.5;
				} else {
					// new cluster
					clusters.push_back(clu);
				}
			}
		}
		////////////////////////////////////////////////////////////////////////

	}

}

#endif	/* MEANSHIFT_H */

