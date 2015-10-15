/*******************************************************************************
 * colorHistMex.cpp - color probability 3d histogram
 *******************************************************************************
 * Add license here...
 *******************************/

#include "mexwrap.h"

typedef pm::MatWrapper FastImage;

class ColorChannel {
public:
	int binCount;
	float min, max;
	float binSize;
	
	ColorChannel() : binCount(1), min(0), max(0), binSize(1){}
	ColorChannel(int N, double m, double M) : binCount(N), min(m), max(M), binSize(range() / binCount){
	}
	
	inline float range() const {
		return max - min;
	}
	
	inline int index(float value) const {
		int idx = std::floor(float(value - min) / binSize);
		return std::max(0, std::min(idx, binCount - 1));
	}
};

/**
 * Usage:
 * 
 * N = colorHistMex( I, options )
 * 
 */
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
	// checking the input
	if (nin > 2) {
		mexErrMsgIdAndTxt("MATLAB:colorhist:maxrhs",
				"Too many arguments, max 2, got %d", nin);
	} else if(nin < 1) {
		mexErrMsgIdAndTxt("MATLAB:colorhist:minrhs",
				"Requires at least an image!");
	}
	
	// read the input image
	const FastImage I(in[0]);
	if(I.channels() != 3){
		mexErrMsgIdAndTxt("MATLAB:colorhist:image", "Unsupported image type, 3 channels only!");
	}
	
	// parameters
	float sigma = 3.0f;
	std::vector<ColorChannel> channels;
	FastImage mask;
	if(nin > 1){
		const mxArray *options = in[1];
		if(!mxIsStruct(options)){
			mexErrMsgIdAndTxt("MATLAB:colorhist:arg", "Second argument should be a struct!");
		}
		mxArray *tmp;
		if(tmp = mxGetField(options, 0, "channels")){
			std::vector<float> data;
			mxLoadVector(data, tmp, "Invalid option 'channels'");
			if(data.size() != 9){
				mexErrMsgIdAndTxt("MATLAB:colorhist:channels", "Requires 3 channels as [bins, min, max]");
			}
			for(int i = 0; i < data.size(); i += 3) {
				channels.push_back(ColorChannel(int(data[i]), data[i + 1], data[i + 2]));
			}
		}
		if(tmp = mxGetField(options, 0, "sigma"))
			sigma = mxCheckedScalar(tmp, "Invalid option 'sigma'");
		if(tmp = mxGetField(options, 0, "mask")){
			mask = FastImage(tmp);
			if(mask.width != I.width || mask.height != I.height || mask.channels() != 1) {
				mexErrMsgIdAndTxt("MATLAB:colorhist:mask", "Invalid mask dimensions!");
			}
		}
	}
	if(channels.empty()){
		channels.push_back(ColorChannel(100, 0.0f, 100.0f)); // L
		channels.push_back(ColorChannel(100, -80.0f, 80.0f)); // a
		channels.push_back(ColorChannel(100, -80.0f, 80.0f)); // b
	}
	const ColorChannel &L = channels[0];
	const ColorChannel &a = channels[1];
	const ColorChannel &b = channels[2];
	
	// the histogram
	FastImage H(L.binCount, a.binCount, mxSINGLE_CLASS, b.binCount);
	if(H.empty()){
		mexErrMsgIdAndTxt("MATLAB:colorhist:memory", "Couldn't allocate space, try with fewer bins!");
	}
	
	// accumulate colors with kernels
	int supportRadius = std::ceil(2 * sigma);
	for(int y = 0; y < I.height; ++y) {
		for(int x = 0; x < I.width; ++x) {
			if(!mask.empty() && !mask.read<bool>(y, x)) continue;
			// read pixel
			pm::Vec3f pixel;
			pm::Vec3i from, to;
			for(int c = 0; c < 3; ++c){
				float v = I.read<float>(y, x, c);
				pixel[c] = v;
				int l = channels[c].index(v);
				from[c] = std::max(0, l - supportRadius);
				to[c] = std::min(l + supportRadius, channels[c].binCount - 1);
			}
			
			// vote for that histogram location (and the 3d neighborhood)
			pm::Vec3f loc = (pixel - pm::Vec3f(L.min, a.min, b.min)).mul(
				pm::Vec3f(1.0f / L.binSize, 1.0f / a.binSize, 1.0f / b.binSize)
			);
			// std::cout << from[0] << ", " << from[1] << ", " << from[2] << "; to=";
			// std::cout << to[0] << ", " << to[1] << ", " << to[2] << "\n";
			for(int m = from[0]; m <= to[0]; ++m) {
				for(int n = from[1]; n <= to[1]; ++n) {
					for(int p = from[2]; p <= to[2]; ++p) {
						pm::Vec3f delta = loc - pm::Vec3f(m, n, p);
						// std::cout << "delta = " << delta[0] << ", " << delta[1] << ", " << delta[2] << "\n";
						float e = std::exp(- delta.dot(delta) / (2.0f * sigma * sigma));
						H.at<float>(m, n, p) += e;
						// std::cout << "at(" << m << ", " << n << ", " << p << ") += " << e << "\n";
					}
				}
			}
			// next pixel
		}
	}
	
	// result
	out[0] = H;
}
