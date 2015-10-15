/*******************************************************************************
 * gb.h - gain + bias implementations
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef GB_H
#define	GB_H

#include <patch.h>
#include <patch/basic.h>
#include <texture.h>

namespace pm {
	template <typename Scalar, int channels, typename Patch>
	inline Vec<Scalar, 3> mean(const Texture *image, const Patch &patch) {
		typedef Vec<Scalar, channels> DataType;
		typedef Vec<Scalar, 3> Vec3;
		Vec3 mu(0.0, 0.0, 0.0);
		const int N = Patch::width();
		for (Point<int> i; i.y < N; ++i.y) {
			for (i.x = 0; i.x < N; ++i.x) {
				DataType p = image->at<DataType>(patch * i);
				mu += Vec3(p[0], p[1], p[2]);
			}
		}
		mu *= 1.0 / Scalar(N * N);
		return mu;
	}

	template <typename Scalar, int channels, typename Patch>
	inline Vec<Scalar, 3> stddev(const Texture *image, const Patch &patch,
			const Vec<Scalar, 3> &mu) {
		typedef Vec<Scalar, channels> DataType;
		typedef Vec<Scalar, 3> Vec3;
		Vec3 sigma(0.0, 0.0, 0.0);
		const int N = Patch::width();
		for (Point<int> i; i.y < N; ++i.y) {
			for (i.x = 0; i.x < N; ++i.x) {
				DataType p = image->at<DataType>(patch * i);
				Vec3 d = Vec3(p[0], p[1], p[2]) - mu;
				sigma += d.mul(d);
			}
		}
		sigma *= 1.0 / Scalar(N * N);
		return Vec3(
				std::sqrt(sigma[0]),
				std::sqrt(sigma[1]),
				std::sqrt(sigma[2])
				);
	}

	template <typename Scalar, int channels>
	inline Vec<Scalar, channels> ranged(const Vec<Scalar, channels> &v,
			const Vec<Scalar, channels> &minv, const Vec<Scalar, channels> &maxv) {
		Vec<Scalar, channels> res;
		for(int i = 0; i < channels; ++i) res[i] = std::min(std::max(v[i], minv[i]), maxv[i]);
		return res;
	}
	
	template <typename Scalar, int channels>
	inline Vec<Scalar, channels> div(const Vec<Scalar, channels> &a, const Vec<Scalar, channels> &b) {
		Vec<Scalar, channels> res;
		for(int i = 0; i < channels; ++i) res[i] = a[i] / b[i];
		return res;
	}
	
	template <typename Scalar = float>
	struct GainBias {
		typedef Scalar Type;
		typedef Vec<Scalar, 3> Vec3;
        typedef unsigned int CacheID;
        
        struct GBData {
            Vec3 mean;
            Vec3 stddev;
            GBData(const Vec3 &m, const Vec3 &s) : mean(m), stddev(s) {}
        };

		static Vec3 minBias;
		static Vec3 maxBias;
		static Vec3 minGain;
		static Vec3 maxGain;
		
		// bias and gain computation
		template <int channels, typename P1, typename P2>
		inline static void compute(const Texture *source, const Texture *target, const P1 &from, const P2 &to,
				Vec3 &gain, Vec3 &bias){
			Vec3 muQ = mean<Scalar, channels>(source, from);
			Vec3 muP = mean<Scalar, channels>(target, to);
			Vec3 siQ = stddev<Scalar, channels>(source, from, muQ);
			Vec3 siP = stddev<Scalar, channels>(target, to, muP);
			gain = ranged(div(siQ, siP), minGain, maxGain);
			bias = ranged(muQ - muP.mul(gain), minBias, maxBias);
		}
        
        template <typename P1, typename P2>
		inline static void fetch(CacheID source, CacheID target, const P1 &from, const P2 &to,
				Vec3 &gain, Vec3 &bias){
            const GBData &gbQ = profile(source, from.y, from.x);
            const GBData &gbP = profile(target, to.y, to.x);
			gain = ranged(div(gbQ.stddev, gbP.stddev), minGain, maxGain);
			bias = ranged(gbQ.mean - gbP.mean.mul(gain), minBias, maxBias);
		}
		
		// bias and gain application
		template <int channels>
		inline static void applyOn(Vec<Scalar, channels> &pixel, const Vec3 &gain, const Vec3 &bias) {
			for (int i = 0; i < 3; ++i) pixel[i] = gain[i] * pixel[i] + bias[i];
		}
        
        static CacheID Target;
        static CacheID Exemplar;
        
        inline static GBData &profile(CacheID id, int y, int x) {
            return caches[id].at(y, x);
        }
        template <int channels>
        inline static void cache(CacheID id, const Texture *img) {
            if(caches.size() <= id) {
                caches.resize(id);
            }
            GBCache &grid = caches[id];
            grid = new GBCache(img->height, img->width);
            typedef BasicPatch<int> Patch;
            for(int y = 0; y < img->height - Patch::width() + 1; ++y) {
                for(int x = 0; x < img->width - Patch::width() + 1; ++x) {
                    Patch p(y, x);
                    GBData &d = grid.at(y, x);
                    d.mean = mean<Scalar, channels>(img, p);
                    d.stddev = stddev<Scalar, channels>(img, p, d.mean);
                }
            }
        }
        
        inline static bool cached() {
            return caches.size() < 2;
        }
        
    private:
        typedef Grid<GBData> GBCache;
        static std::vector<GBCache> caches;
	};
	
	#define DEFINE_PROP(name, value) \
	template <typename Scalar> \
	Vec<Scalar, 3> GainBias<Scalar>::name(value, value, value)
	// static stuff:
	DEFINE_PROP(minBias, 0);
	DEFINE_PROP(maxBias, 0);
	DEFINE_PROP(minGain, 1);
	DEFINE_PROP(maxGain, 1);
   
    template <typename Scalar>
    std::vector<typename GainBias<Scalar>::GBCache> GainBias<Scalar>::caches(2);
    template <typename Scalar>
    typename GainBias<Scalar>::CacheID GainBias<Scalar>::Target(0);
    template <typename Scalar>
    typename GainBias<Scalar>::CacheID GainBias<Scalar>::Exemplar(1);
}

#endif	/* GB_H */

