/*******************************************************************************
 * jigsawMex.cpp - smart init implementation for matlab using partial xcorr
 *******************************************************************************
 * Add license here...
 *******************************/

#include "algebra.h"
#include "mexwrap.h"
#include "rng.h"
#include <vector>

using namespace pm;

//typedef MatWrapper MatXD;
typedef unsigned int uint;
//typedef Point<int> Point2i;

struct Block {
	Point2i from;
	Point2i to;

	explicit Block(const Point2i &from = Point2i(), const Point2i &to = Point2i()) : from(from), to(to) {
	}
	
	inline void shift(const Point2i delta) {
		from = from + delta;
		to = to + delta;
	}

	inline int width() const {
		return to.x - from.x + 1;
	}

	inline int height() const {
		return to.y - from.y + 1;
	}

	inline Point2i size() const {
		return Point2i(width(), height());
	}

	inline bool contains(const Point2i &p) const {
		return p.x >= from.x && p.y >= from.y && p.x <= to.x && p.y <= to.y;
	}
	
	inline Block expand(const Point2i &margin) const {
		return Block(from - margin, to + margin);
	}
};

struct JigsawBlockIterator {
	const Point2i fullSize;
	const Point2i blockSize;
	
	JigsawBlockIterator(const Point2i &full, const Point2i &blk) : fullSize(full), blockSize(blk) {
		i = 0;
		numBlocks = Point2i(
			std::ceil(float(fullSize.x) / blockSize.x),
			std::ceil(float(fullSize.y) / blockSize.y)
		);
	}
	
	inline Point2i idx() const {
		return Point2i(i % numBlocks.x, i / numBlocks.y);
	}
	
	inline Block fullBlock() const {
		Point2i pos = idx();
		Point2i from = pos.mult(blockSize);
		Point2i to = from + blockSize - Point2i(1, 1);
		return Block(from, to);
	}
	
	inline Block block() const {
		Point2i pos = idx();
		Point2i from = pos.mult(blockSize);
		Point2i to(
			std::min(fullSize.x, from.x + blockSize.x) - 1, 
			std::min(fullSize.y, from.y + blockSize.y) - 1
		);
		return Block(from, to);
	}
	
	inline Block operator *() const {
		return block();
	}
	
	inline operator bool() const {
		return i < numBlocks.area();
	}
	
	inline JigsawBlockIterator &operator ++() {
		++i;
		return *this;
	}
	
	inline bool isFirst() const {
		return i == 0;
	}
	
private:
	int i;
	Point2i numBlocks;
};

struct ImageBlockView {
	const MatXD *image;
	const Block block;
	const MatXD *mask;
	ImageBlockView(const MatXD *img, const Block &blk, const MatXD *m = NULL) : image(img), block(blk), mask(m) {
		// allocate the mean
		mean.resize(image->channels(), 0);
		// compute the mean
		for(int ch = 0, num_ch = img->channels(); ch < num_ch; ++ch) {
			int N = 0;
			for(int x = blk.from.x; x <= blk.to.x; ++x) {
				for(int y = blk.from.y; y < blk.to.y; ++y) {
					// mean within the block, within the image, only for the valid mask parts
					if(contains(y, x)){
						mean[ch] += image->read<float>(y, x, ch);
						++N;
					}
				}
			}
			mean[ch] /= N;
		}
	}
	
	inline float get(int dy, int dx, int ch) const {
		return image->read<float>(block.from.y + dy, block.from.x + dx, ch) - mean[ch];
	}
	
	inline bool contains(int y, int x) const {
		return y >= 0 && y < image->height && x >= 0 && x < image->width
				&& (mask == NULL || mask->read<bool>(y, x));
	}
	
	inline bool contains(const Point2i &delta) const {
		return contains(block.from.y + delta.y, block.from.x + delta.x);
	}
	
	inline float meanValue(int ch) const {
		return mean[ch];
	}
	
private:
	std::vector<float> mean;
};

float correlation(const ImageBlockView &v1, const ImageBlockView &v2, const Point2i &delta) {
	float sum = 0.0f;
	for(int ch = 0, num_ch = v1.image->channels(); ch < num_ch; ++ch) {
		const Point2i size = v1.block.size();
		for(int dx = 0; dx < size.x; ++dx) {
			for(int dy = 0; dy < size.y; ++dy) {
				Point2i d1 = Point2i(dx, dy);
				Point2i d2 = d1 + delta;
				// mean within the block, within the image, only for the valid mask parts
				if(v1.contains(d1) && v2.contains(d2)){
					sum += v1.get(d1.y, d1.x, ch) * v2.get(d2.y, d2.x, ch);
				}
			}
		}
	}
	return sum;
}

/**
 * Usage:
 * 
 * [ T ] = jigsawmex( T, S, G, seed )
 */
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
	// checking the input
	if (nin != 4) {
		mexErrMsgIdAndTxt("MATLAB:jigsaw:rhs",
				"Requires 4 arguments, got //d", nin);
	}

	// checking the output
	if (nout != 1) {
		mexErrMsgIdAndTxt("MATLAB:jigsaw:lhs", "I have 1 only output value!");
	}

	// arguments
	const MatXD T0(in[0]);
	const MatXD S(in[1]);
	const MatXD G(in[2]);
	unsigned int s = mxCheckedScalar(in[3], "Invalid seed parameter!");

	// seed the random part
	seed(s);
	std::cout << "Seeded with " << s << "\n";

	// block size P
	Point2i g1(G.read<int>(0, 1), G.read<int>(0, 0));
	Point2i g2(G.read<int>(1, 1), G.read<int>(1, 0));
	Point2i P = g1.area() >= g2.area() ? g1 : g2;

	// blocks position
	int srcDimY = std::floor(float(S.height) / P.y);
	int srcDimX = std::floor(float(S.width) / P.x);
	
	// result data
	MatXD T(T0.height, T0.width, T0.type()); // the result data
	MatXD mask(T0.height, T0.width, mxLOGICAL_CLASS); // the processing mask
	
	std::cout << "Starting.\n";

	// build the puzzle, row by row, block by block
	JigsawBlockIterator it(Point2i(T.width, T.height), P);
	for(; it; ++it){
		// target block
		const Block trgBlock = *it;
		
		std::cout << "Block from " << trgBlock.from.y << "/" << trgBlock.from.x;
		std::cout << " to " << trgBlock.to.y << "/" << trgBlock.to.x << ":\n";
		
		// source block
		int x = uniform(unif01, 0, srcDimX - 1);
		int y = uniform(unif01, 0, srcDimY - 1);
		Block srcBlock;
		{
			Point2i from(x * P.x, y * P.y);
			Point2i to(std::min(from.x + P.x, S.width), std::min(from.y + P.y, S.height));
			srcBlock = Block(from, to);
		}
		
		// if it's the first block, just copy it since no alignment is possible
		if(!it.isFirst()) {
			// the block views
			ImageBlockView trgView(&T, it.fullBlock().expand(P), &mask);
			Block srcFullBlock(Point2i(x - 1, y - 1).mult(P), Point2i(x + 2, y + 2).mult(P) - Point2i(1, 1));
			ImageBlockView srcView(&S, srcFullBlock);
			
			// get the best reasonable shift
			Point2i from = Point2i::max(Point2i(0, 0), srcBlock.from - P);
			Point2i to = Point2i::min(Point2i(S.width, S.height), srcBlock.to + P) - Point2i(1, 1);
			Block shifts(from, to);
			float bestCorr = -1.0f;
			Point2i bestShift;
			for(int dy = shifts.from.y; dy <= shifts.to.y; ++dy) {
				for(int dx = shifts.from.x; dx <= shifts.to.x; ++dx) {
					Point2i shift = srcBlock.from - Point2i(dx, dy);
					float corr = correlation(trgView, srcView, shift);
					std::cout << "d(" << shift.y << "/" << shift.x << ") -> corr=" << corr << "\n";
					if(corr > bestCorr) {
						bestCorr = corr;
						bestShift = shift;
					}
				}
			}
			std::cout << "Best shift = " << bestShift.y << "/" << bestShift.x << " for corr=" << bestCorr << "\n";
			
			// we use the best shift
			srcBlock.shift(bestShift);
		}
		
		// copy the block from the source to the result
		// do that using locality (help the cache, please!)
		for(int ch = 0; ch < T.channels(); ++ch) {
			for(int tX = trgBlock.from.x, sX = srcBlock.from.x; tX <= trgBlock.to.x; ++tX, ++sX ){
				for(int tY = trgBlock.from.y, sY = srcBlock.from.y; tY <= trgBlock.to.y; ++tY, ++sY ){
					T.update(tY, tX, ch, S.read<float>(sY, sX, ch));
					mask.update(tY, tX, true);
				}
			}
		}
	}
	
	out[0] = T;
}
