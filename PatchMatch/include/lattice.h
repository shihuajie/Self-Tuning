/*******************************************************************************
 * lattice.h - lattice extraction implementation
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef LATTICE_H
#define	LATTICE_H

#define _USE_MATH_DEFINES

#include "algebra.h"
#include <algorithm>
#include <map>
#include <vector>

namespace pm {
	
	inline bool operator==(const Point2i &p1, const Point2i &p2) {
		return p1.x == p2.x && p1.y == p2.y;
	}
	struct point_order {
		bool operator()(const pm::Point2i &p1, const pm::Point2i &p2) {
			return p1.x == p2.x ? p1.y < p2.y : p1.x < p2.x;
		}
	};
	
	struct Bounds {
		Point2i min;
		Point2i max;
		explicit Bounds(const Point2i &m1 = Point2i(), const Point2i &m2 = Point2i()) : min(m1), max(m2) {}
		inline int width() const {
			return max.x - min.x + 1;
		}
		inline int height() const {
			return max.y - min.y + 1;
		}
		inline int size() const {
			return width() * height();
		}
		inline bool contains(const Point2i &p) const {
			return p.x >= min.x && p.y >= min.y && p.x <= max.x && p.y <= max.y;
		}
		
		class MultipleIterator {
		public:
			MultipleIterator(const Bounds *p, const Point2i &d, const Point2i &s) : parent(p), delta(d), shift(s), i(0) {
				if(!isValid()) {
					selectNext();
				}
			}
			
			inline bool isValid() const {
				return parent->contains(shift + delta * i);
			}
			inline operator bool() const {
				return isValid();
			}
			inline Point2i operator*() const {
				return shift + delta * i;
			}
			inline MultipleIterator &operator++() {
				// pre-increment
				selectNext();
				return *this;
			}
			inline MultipleIterator operator++(int) {
				// post-increment
				MultipleIterator it = *this; // copy current state
				selectNext(); // this changes
				return it; // we return the previous state
			}
			
		protected:
			inline void selectNext() {
				if(i >= 0) {
					++i;
					if(!isValid()){
						i = -1;
					}
				} else {
					--i;
				}
			}
			
		private:
			const Bounds *parent;
			Point2i delta;
			Point2i shift;
			int i;
		};
		
		inline MultipleIterator multiple_begin(const Point2i &delta, const Point2i &shift = Point2i()) const {
			return MultipleIterator(this, delta, shift);
		}
	};
	
	struct DistanceRange {
		float min, max;
		explicit DistanceRange(float from = 5.0f, float to = FLT_MAX) : min(from), max(to) {}
		
		template <typename D>
		bool contains(Point<D> p) const {
			float dist = float(sqrt((double)p.dot(p)));
			return dist >= min && dist <= max;
		}
	};
	
	struct EnergyTerm {
		float sum;
		unsigned int N;
		Point2i shift;
		
		explicit EnergyTerm(float s = 0, unsigned int n = 0, Point2i sh = Point2i()) : sum(s), N(n), shift(sh){}
		
		float value() const {
			if(N == 0) return 0.0f;
			float dist = sqrt((double)shift.dot(shift));
			return sum / (N * dist);
		}
		operator float() const {
			return value();
		}
	};
	
	struct Generator {
		Point2i point;
		float energy;
		explicit Generator(const Point2i &p = Point2i(), float e = 0.0f) : point(p), energy(e) {}
	};
	
	struct GeneratorPair {
		Point2i p1, p2;
		float e1, e2;
		explicit GeneratorPair(const Point2i &p1 = Point2i(), const Point2i &p2 = Point2i(), 
			float e1 = 0.0f, float e2 = -1.0f) : p1(p1), p2(p2), e1(e1), e2(e2 < 0 ? e1 : e2) {}
	};
	
	bool compare_gen(const Generator &g1, const Generator &g2) {
		return g1.energy > g2.energy;
	}
	
	inline float angleDiff(const Point2i &p1, const Point2i &p2) {
		const Point2d d1(p1), d2(p2); // work with double values to be safe
		double sqNorms = d1.sqNorm() * d2.sqNorm();
		double cosDelta = d1.dot(d2) / std::sqrt(sqNorms);
		
		// angle errors may lead to nan
		if(std::abs(cosDelta) >= 1.0) return 0.0f;
		double delta = std::acos(cosDelta);
		if(_isnan(delta)) return 0.0f;
		
		// take only the positive <90° difference
		delta = std::abs(delta); // use only positive angles
		delta = delta > M_PI_2 ? M_PI - delta : delta; // get the <= 90° angle diff
		return float(delta);
	}
	
	class Lattice {
	public:
		
		enum Measure {
			Precision,
			Recall,
			FMeasure,
			GMeasure,
			Value
		};
		
		Lattice(const std::vector<Point2f> &offsets, DistanceRange r, 
		float sigma = 3.0f, float penFactor = 1.0f)
				: distRange(r), map(NULL) {
			if(offsets.empty()) return; // we cannot do anything
			
			// find the bounds
			Point2f first = offsets[0];
			Point2i min = first.floor(), max = first.ceil();
			for(unsigned int i = 1, n = offsets.size(); i < n; ++i) {
				min = Point2i::min(min, offsets[i].floor());
				max = Point2i::max(max, offsets[i].ceil());
			}
			int supportRadius = std::ceil(2 * sigma);
			Point2i support(supportRadius, supportRadius);
			bounds = Bounds(min - support, max + support);
			
			// use min-shift and generate the map
			map = new float[bounds.size()](); //< init to 0.0
			for(unsigned int i = 0, n = offsets.size(); i < n; ++i) {
				const Point2f &p = offsets[i];
				// add the support values
				Point2i from = p.floor() - support;
				Point2i to = p.ceil() + support;
				for(int y = from.y; y <= to.y; ++y) {
					for(int x = from.x; x <= to.x; ++x) {
						Point2f delta = p - Point2f(x, y);
						at(y, x) += std::exp(- delta.dot(delta) / (2.0f * sigma));
					}
				}
			}
			
			// computation of the mean non-zero energy
			float nonZeroSum = 0.0f;
			unsigned int nnz = 0;
			for(int y = bounds.min.y; y <= bounds.max.y; ++y) {
				for(int x = bounds.min.x; x <= bounds.max.x; ++x) {
					float value = at(Point2f(x, y));
					if(value > 0){
						nonZeroSum += value;
						++nnz;
					}
				}
			}
			penalty = nonZeroSum / nnz * penFactor;
			total = nonZeroSum;
		}
		
		virtual ~Lattice() {
			if (map != NULL) delete[] map;
		}
		
		inline const float &at(const Point2i &p) const {
			return at(p.y, p.x);
		}
		inline const float &at(int y, int x) const {
			return map[bounds.width() * (y - bounds.min.y) + (x - bounds.min.x)];
		}
		
		// the energy function to minimize to get the lattice generators
		inline EnergyTerm energy(const Point2i &x1, const Point2i &x2 = Point2i()) const {
			float sum = 0.0f;
			unsigned int N = 0;
			if(x2.x == 0 && x2.y == 0) {
				Bounds::MultipleIterator it = bounds.multiple_begin(x1);
				while(it){
					// get the multiple (non-zero)
					Point2i d = *it;
					// except for the origin which is singular
					if(d.x != 0 || d.y != 0){
						float c = at(d); // add its contribution
						// non-null or we get a penalty!
						if(c > 0.0f){
							sum += c;
						} else {
							sum -= penalty;
						}
						++N; // increase the valid counts
					}
					++it; // next multiple
				}
				return EnergyTerm(sum, N, x1);
				
			} else {
				// general 2d generator case
				// first multiple iterator
				Bounds::MultipleIterator it1 = bounds.multiple_begin(x1);
				while(it1){
					const Point2i d = *it1;
					
					// second multiple iterator
					Bounds::MultipleIterator it2 = bounds.multiple_begin(x2, d);
					while(it2) {
						Point2i z = *it2;
						if(z.x != 0 || z.y != 0) {
							float c = at(z); // add its contribution
							// non-null or we get a penalty!
							if(c > 0.0f){
								sum += c;
							} else {
								sum -= penalty;
							}
							++N; // increase the valid counts
						}
						++it2;
					}
					
					++it1; // next multiple
				}
				return EnergyTerm(sum, N, x2);
			}
		}
		
		// precision = #positive coverage / #pos+neg coverage
		inline float precision(const EnergyTerm &e) const {
			if(e.N > 0 && e.sum > 0)
				return e.sum / e.N;
			return 0.0f;
		}
		
		// recall = #coverage / max#coverage
		inline float recall(const EnergyTerm &e) const {
			return std::max(e.sum, 0.0f) / total;
		}
		
		// F-beta = (1+beta²) precision * recall / (beta² precision + recall)
		inline float Fmeasure(const EnergyTerm &e, float beta = 1.0f) const {
			float betaSq = beta * beta;
			float P = precision(e);
			float R = recall(e);
			if(P < 1e-8 || R < 1e-8) return 0.0f;
			return (1 + betaSq) * P * R / (betaSq * P + R); 
		}
		
		// G-measure = sqrt(precision * recall)
		inline float Gmeasure(const EnergyTerm &e) const {
			return precision(e) * recall(e);
		}
		
		inline float measure(const EnergyTerm &e, Measure m = FMeasure, float beta = 0.1f) const {
			switch(m){
				case Precision: return precision(e);
				case Recall:	return recall(e);
				case FMeasure:	return Fmeasure(e, beta);
				case GMeasure:	return Gmeasure(e);
				case Value:		return e.value();
				default:
					std::cerr << "Invalid measure: " << m << "\n";
					return 0.0f;
			}
		}
		
		// since g1 and -g1 have the same energy, we only consider the right quadrants
		inline Generator getFirstGenerator(Measure m = FMeasure, float beta = 0.1f) const {
			float maxEnergy = 0;
			Point2i maxOffset;
			for(int y = bounds.min.y; y <= bounds.max.y; ++y) {
				for(int x = 0; x <= bounds.max.x; ++x) {
					Point2i d(x, y);
					// we skip the ones out of range
					if(!distRange.contains(d)) continue;
					// we get the energy
					EnergyTerm E = energy(d);
					float e = measure(E, m, beta);
					if(e > maxEnergy) {
						maxEnergy = e;
						maxOffset = d;
					}
				}
			}
			return Generator(maxOffset, maxEnergy);
		}
		inline std::vector<Generator> getNFirstGenerators(unsigned int N = 10, 
		Measure m = FMeasure, float beta = 0.1f) const {
			std::vector<Generator> gen;
			for(int y = bounds.min.y; y <= bounds.max.y; ++y) {
				for(int x = 0; x <= bounds.max.x; ++x) {
					Point2i d(x, y);
					// we skip the ones out of range
					if(!distRange.contains(d)) continue;
					// we get the energy
					EnergyTerm E = energy(d);
					float e = measure(E, m, beta);
					if(e > 0) {
						gen.push_back(Generator(d, e));
					}
				}
			}
			// sort by energy
			std::sort(gen.begin(), gen.end(), compare_gen);
			// keep only N
			if(gen.size() > N) gen.resize(N);
			return gen;
		}
		
		// since g2 and -g2 have the same energy, we only consider the top quadrants
		inline Generator getSecondGenerator(const Point2i &g1, float minAngle = M_PI / 16,
			Measure m = FMeasure, float beta = 0.1f) const {
			float maxEnergy = 0;
			Point2i maxOffset;
			for(int y = 0; y <= bounds.max.y; ++y) {
				for(int x = bounds.min.x; x <= bounds.max.x; ++x) {
					Point2i d(x, y);
					// we skip the ones out of range
					if(!distRange.contains(d)) continue;
					// we also skip the offsets which are almost multiples of g1
					if(angleDiff(g1, d) < minAngle) continue; // too close from a multiple of g1
					// we get the energy
					EnergyTerm et = energy(g1, d);
					float e = measure(et, m, beta);
					if(e > maxEnergy) {
						maxEnergy = e;
						maxOffset = d;
					}
				}
			}
			return Generator(maxOffset, maxEnergy);
		}
		
		// get g1 and g2 by maximizing an energy together (deadly slow)
		inline GeneratorPair getBestGeneratorPair(float minAngle = M_PI / 16,
			Measure measure = FMeasure, float beta = 0.1f) const {
			float maxEnergy = -1.0f;
			GeneratorPair maxPair;
			for(int y = bounds.min.y; y <= bounds.max.y; ++y) {
				for(int x = 0; x <= bounds.max.x; ++x) {
					Point2i d(x, y);
					// we skip the ones out of range
					if(!distRange.contains(d)) continue;
					// std::cout << "@(" << y << ", " << x << "): ";
					// we get the energy
					Generator g2 = getSecondGenerator(d, minAngle, measure, beta);
					// std::cout << g2.point.y << ", " << g2.point.x << " => " << g2.energy << "\n";
					if(g2.energy > maxEnergy) {
						maxEnergy = g2.energy;
						maxPair = GeneratorPair(d, g2.point, g2.energy);
					}
				}
			}
			return maxPair;
		}
		
		// get g1 and g2 by maximizing for g1 and keeping the N best
		inline GeneratorPair approxBestGeneratorPair(unsigned int N = 100,
			float minAngle = M_PI / 16, Measure measure = FMeasure, float beta = 0.1f) const {
			typedef std::vector<Generator> GenList;
			GenList bestG1 = getNFirstGenerators(N, measure, beta);
			float maxEnergy = -1.0f;
			GeneratorPair maxPair;
			for(GenList::iterator it = bestG1.begin(); it != bestG1.end(); ++it) {
				Generator g1 = *it;
				// we get the energy
				Generator g2 = getSecondGenerator(g1.point, minAngle, measure, beta);
				if(g2.energy > maxEnergy) {
					maxEnergy = g2.energy;
					maxPair = GeneratorPair(g1.point, g2.point, g1.energy, g2.energy);
				}
			}
			return maxPair;
		}
		
		inline Generator getStableGenerator(const Point2i &g1, 
				float minAngle = M_PI / 16, Measure m = FMeasure, float beta = 0.1f, unsigned int n = 10) const {
			std::map<pm::Point2i, bool, point_order> seeds;
			Generator seed; 
			seeds[g1] = true;
			for(unsigned int i = 0, N = 2 * n; i < N; ++i) {
				Generator gen = getSecondGenerator(seed.point, minAngle, m, beta);
				seed = gen;
				// if we've already seen it, we're done
				if(seeds[gen.point]){
					return seed;
				} else {
					seeds[gen.point] = true;
				}
			}
			return seed;
		}
		
		inline const Bounds &boundaries() const {
			return bounds;
		}
		
	protected:
		inline float &at(int y, int x) {
			return map[bounds.width() * (y - bounds.min.y) + (x - bounds.min.x)];
		}
				
	private:
		DistanceRange distRange;
		Bounds bounds;
		float *map;
		float penalty, total;
	};
	
}

#endif	/* LATTICE_H */

