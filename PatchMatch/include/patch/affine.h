/*******************************************************************************
 * affine.h - affine patch implementation
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef AFFINE_H
#define	AFFINE_H

#include "basic.h"

// we declare sincosf to help IDEs
void sincosf(float, float*, float*);

#ifndef _GNU_SOURCE 
inline void sincosf(float a, float* s, float* c) {
	*s = std::sin(a);
	*c = std::cos(a);
}
#endif

namespace pm {

	/**
	 * \brief Type of mirroring (scale inversion)
	 */
	enum Mirroring {
		None = 0,
		MirrorX = 1,
		MirrorY = 2,
		MirrorXY = 3
	};

	/**
	 * \brief Properties of the space of affine patches
	 */
	template <int AngleSteps = 6 >
	struct AffineProperties {
		static const int angleSteps = AngleSteps;
		static const int angleCount = 2 * angleSteps + 1;

		// Parameters
		float minScale;
		float maxScale;
		float maxAngle;
		Mirroring mirroring;
		bool homogeneousScaling;

		/**
		 * \brief Default constructor
		 */
		AffineProperties() {
			minScale = 0.9f;
			maxScale = 1.2f;
			maxAngle = M_PI_4;
			mirroring = MirrorX;
			homogeneousScaling = false;
			// generate tables
			initTables();
		}

		// The tables (declarations)
		float angleTable[angleCount];
		float cosTable[angleCount];
		float sinTable[angleCount];

		/**
		 * \brief To be called when maxAngle is changed
		 */
		void initTables() {
			// init angles + cos/sin
			angleTable[AngleSteps] = 0.0;
			cosTable[AngleSteps] = 1.0;
			sinTable[AngleSteps] = 0.0;
			for (int i = 1; i <= AngleSteps; ++i) {
				float alpha = maxAngle / AngleSteps * i;
				angleTable[AngleSteps + i] = alpha;
				angleTable[AngleSteps - i] = -alpha;
				float c, s;
				sincosf(alpha, &s, &c);
				sinTable[AngleSteps + i] = s;
				sinTable[AngleSteps - i] = -s;
				cosTable[AngleSteps + i] = cosTable[AngleSteps - i] = c;
			}
		}
	};

	/**
	 * \brief The default properties for affine patches
	 */
	static AffineProperties<> DefaultProperties;

	/**
	 * \brief Affine patch representation
	 */
	template <int PatchWidth = DEFAULT_PATCH_SIZE>
	class AffinePatch : public BasicGrid < AffinePatch<PatchWidth> > {
	public:

		typedef int TableIndex;
		typedef BasicPatch<int, PatchWidth> OriginalPatchType;
		typedef Point<float> PixLoc;
		typedef Point<int> Index;
		typedef float Coherence;
		// dimensions: x,y,a,sx (, sy if homogeneous + no mirroring) 

		inline static int dimensions() {
			return props.homogeneousScaling && props.mirroring == None ? 4 : 5;
		}

		/// The Patch size

		inline static int width(int newSize = 0) {
			if (newSize != 0) throw "Cannot change static width!";
			return PatchWidth;
		}

		inline static float radius() {
			return 0.5f * width();
		}

		// The Patch Space properties
		static AffineProperties<> props;
		
		// The patch coordinates
		float x, y;
		// The patch transformation
		TableIndex a;
		float sx, sy;
		int angleIndex;
		bool operator==(const AffinePatch<PatchWidth> &p) const {
			return x == p.x && y == p.y && a == p.a && sx == p.sx && sy == p.sy;
		}

		/**
		 * The transformation within the patch center system is:
		 * 
		 * (x', y') = T * (x, y)
		 * 
		 * where T = [cos -sin  x [sx 0  = [sx cos  -sy sin
		 *			  sin  cos]    0 sy]    sx sin   sy cos]
		 */

		inline static PixLoc pxloc(float x, float y, float sx, float sy, int a, int dx, int dy) {
			// translation to the center
			float px = dx; // - radius();
			float py = dy; // - radius();
			// transformation paramaters
			float cosA = props.cosTable[a];
			float sinA = props.sinTable[a];
			// actual transformation
			// 1 = scale / mirror
			px *= sx;
			py *= sy;
			// 2 = rotate
			float tx = cosA * px - sinA * py;
			float ty = sinA * px + cosA * py;
			// 3 = translate
			return PixLoc(tx + x, ty + y);
		}

		bool withinFrame(const Image *frame) const {
			const int maxX = frame->cols;
			const int maxY = frame->rows;
			// basics
			if (x < 0 || x > maxX - width() || y < 0 || y > maxY - width()) return false;
			// top left
			PixLoc px = pxloc(x, y, sx, sy, a, 0, 0);
			if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;
			// top right
			px = pxloc(x, y, sx, sy, a, 0, width() - 1);
			if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;
			// top left
			px = pxloc(x, y, sx, sy, a, width() - 1, 0);
			if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;
			// top left
			px = pxloc(x, y, sx, sy, a, width() - 1, width() - 1);
			if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;

			// everything is ok!
			return true;
		}

		inline void points(PixLoc pts[4]) const {
			pts[0] = pxloc(x, y, sx, sy, a, 0, 0);
			pts[1] = pxloc(x, y, sx, sy, a, 0, width() - 1);
			pts[2] = pxloc(x, y, sx, sy, a, width() - 1, width() - 1);
			pts[3] = pxloc(x, y, sx, sy, a, width() - 1, 0);
		}

		inline float cx() const {
			return x + radius();
		}

		inline float cy() const {
			return y + radius();
		}

		void debug() const {
			std::cout << "x=" << x << ", y=" << y << ", sx=" << sx << ", sy=" << sy << "\n";
			std::cout << "\ta=" << a << " => cos=" << props.cosTable[a] << ", sin=" << props.sinTable[a] << "\n";
		}

		template <typename Storage>
		inline void store(Storage &out, int channels) const {
			switch (channels) {
				case 6:
					out[5] = width();
				case 5:
					out[4] = sy;
				case 4:
					out[3] = sx;
				case 3:
					out[2] = props.angleTable[a];
				case 2:
					out[1] = x;
					out[0] = y;
					break;
				default:
					std::cerr << "Unknown storage of " << channels;
					std::cerr << " channels!\n";
			}
		}

		inline void load(const float *source, int channels, int offset) {
			if(channels != dimensions()){
				std::cerr << "Warning: loading an incomplete source of " << channels;
				std::cerr << "channels (current=" << dimensions() << " channels)\n";
			}
			switch (channels) {
				case 5:
					sy = source[offset * 4];
				case 4:
					sx = source[offset * 3];
				case 3:
				{
					float ang = source[offset * 2];
					int best = 0;
					float bestDiff = 10.0f;
					for (int k = 0; k < props.angleCount; ++k) {
						float diff = std::abs(ang - props.angleTable[k]);
						if (diff < bestDiff) {
							bestDiff = diff;
							best = k;
						}
					}
					a = best;
				}
				case 2:
					y = source[0];
					x = source[offset];
					break;
				default:
					std::cerr << "Unknown source of " << channels;
					std::cerr << " channels!\n";
			}
		}
		
		inline PixLoc transform(int py, int px) const {
			return pxloc(x, y, a, sx, sy, px, py);
		}
		inline PixLoc transform(const Index &i) const {
			return transform(i.y, i.x);
		}
		inline PixLoc operator *(const Index &i) const {
			return transform(i);
		}
	};

	template <int PatchSize>
	AffineProperties<> AffinePatch<PatchSize>::props = DefaultProperties;

	// implementation for patches of dynamic size

	template <>
	int AffinePatch<0>::width(int newSize) {
		return BasicPatch<int, 0>::width(newSize); // delegate, to be smart
	}

	/**
	 * \brief Dynamic sized affine patch
	 */
	typedef AffinePatch<0> AffinePatchX;

#define TOO_MANY_ITERATIONS 1000

	template <int PatchWidth>
	inline void randomInit(RNG rand, const Image *parent,
			AffinePatch<PatchWidth> &patch) {
		typedef AffinePatch<PatchWidth> Patch;

		float maxX = parent->cols - Patch::width();
		float maxY = parent->rows - Patch::width();
		int it = 0;
		do {
			if (it++ > TOO_MANY_ITERATIONS) {
				std::cerr << "Spent " << it << " iterations to initialize patch (width=" << Patch::width;
				std::cerr << " in [" << parent->cols << "x" << parent->rows << "])\n";
				patch.sx = 1.0f;
				patch.sy = 1.0f;
				patch.a = Patch::props.angleSteps;
				break;
			}
			patch.x = uniform(rand, 0.0f, maxX);
			patch.y = uniform(rand, 0.0f, maxY);
			patch.sx = uniform(rand, Patch::props.minScale, Patch::props.maxScale);
			if (Patch::props.homogeneousScaling)
				patch.sy = patch.sx;
			else
				patch.sy = uniform(rand, Patch::props.minScale, Patch::props.maxScale);
			patch.a = uniform(rand, 0, 2 * Patch::props.angleSteps);
			switch (Patch::props.mirroring) {
				case MirrorX:
					if (bernoulli(rand)) patch.sx *= -1;
					break;
				case MirrorY:
					if (bernoulli(rand)) patch.sy *= -1;
					break;
				case MirrorXY:
					if (bernoulli(rand)) patch.sx *= -1;
					if (bernoulli(rand)) patch.sy *= -1;
					break;
				default:
					break;
			}
		} while (!patch.withinFrame(parent));
	}

	template <int PatchWidth>
	inline bool random(RNG rand, const Image *parent,
			const AffinePatch<PatchWidth> &oldPatch,
			AffinePatch<PatchWidth> &newPatch,
			int windowSize) {
		typedef AffinePatch<PatchWidth> Patch;

		newPatch.x = uniform(rand,
				std::max(0.0f, oldPatch.x - windowSize),
				std::min(float(parent->cols - Patch::width()), oldPatch.x + windowSize)
				);
		newPatch.y = uniform(rand,
				std::max(0.0f, oldPatch.y - windowSize),
				std::min(float(parent->rows - Patch::width()), oldPatch.y + windowSize)
				);
		newPatch.sx = uniform(rand, Patch::props.minScale, Patch::props.maxScale);
		if (Patch::props.homogeneousScaling)
			newPatch.sy = newPatch.sx;
		else
			newPatch.sy = uniform(rand, Patch::props.minScale, Patch::props.maxScale);
		newPatch.a = uniform(rand, 0, 2 * Patch::props.angleSteps);
		switch (Patch::props.mirroring) {
			case MirrorX:
				if (bernoulli(rand)) newPatch.sx *= -1;
				break;
			case MirrorY:
				if (bernoulli(rand)) newPatch.sy *= -1;
				break;
			case MirrorXY:
				if (bernoulli(rand)) newPatch.sx *= -1;
				if (bernoulli(rand)) newPatch.sy *= -1;
				break;
			default:
				break;
		}
		// the random patch may not be valid!
		return newPatch.withinFrame(parent);
	}
    
    template <int PatchWidth>
	inline bool aligned(RNG rand, const Image *parent,
			const AffinePatch<PatchWidth> &oldPatch,
            AffinePatch<PatchWidth> &newPatch, 
            const Point2f &g1, const Point2f &g2, float jitter) {
        typedef AffinePatch<PatchWidth> Patch;
		typedef BasicPatch<float, PatchWidth> SimplePatch;
		SimplePatch::width(Patch::width());
		
		// pretend to be basic, no problem anymore ...
		SimplePatch tmpPatch(oldPatch.x, oldPatch.y);
		SimplePatch tmpNewPatch;
		using patch::aligned;
		aligned(rand, parent, tmpPatch, tmpNewPatch, g1, g2, jitter);
		
		// transfer to new patch
		newPatch = oldPatch;
		newPatch.x = tmpNewPatch.x;
		newPatch.y = tmpNewPatch.y;
		
		// maybe we're unlucky ...
        return newPatch.withinFrame(parent);
    }

	template <int PatchWidth>
	inline void deltaPatch(
			const AffinePatch<PatchWidth> &patch,
			AffinePatch<PatchWidth> &delta,
			int dy, int dx) {
		typedef AffinePatch<PatchWidth> Patch;
		typedef typename Patch::PixLoc PixLoc;
		// copy everything
		delta = patch;
		// update position
		PixLoc dp = Patch::pxloc(patch.x, patch.y, patch.sx, patch.sy, patch.a, dx, dy);
		delta.x = dp.x;
		delta.y = dp.y;
	}

	template <int PatchWidth>
	inline bool isWithin(const Image *parent, const AffinePatch<PatchWidth> &patch) {
		// typedef typename AffinePatch<D, PatchWidth, AngleSteps> Patch;
		return patch.withinFrame(parent);
		// return patch.x >= 0 && patch.y >= 0 && patch.x < maxX && patch.y < maxY;
	}

	/// Coherence: max(0, 1 - dist(c', c+d) )

	template <int PatchWidth>
	inline typename AffinePatch<PatchWidth>::Coherence coherence(const AffinePatch<PatchWidth> &p1,
			const AffinePatch<PatchWidth> &p2, int dy, int dx) {
		typedef AffinePatch<PatchWidth> Patch;
		typedef typename Patch::PixLoc PixLoc;
		int radius = int(Patch::radius());
		PixLoc c1pd = Patch::pxloc(p1.x, p1.y, p1.sx, p1.sy, p1.a, radius + dx, radius + dy);
		PixLoc c2 = Patch::pxloc(p2.x, p2.y, p2.sx, p2.sy, p2.a, radius, radius);
		PixLoc cdiff = c2 - c1pd;
		return std::max(0.0f, 1.0f - std::sqrt(cdiff.dot(cdiff)));
	}
}

#endif	/* AFFINE_H */

