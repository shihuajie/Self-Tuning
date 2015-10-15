/*******************************************************************************
 * voteMex.cpp - voting matlab function implementation
 *******************************************************************************
 * Add license here...
 *******************************/

// other stuff
#include <mexutil.h>
// our stuff
#include <distance.h>
#include <filter.h>
#include <nnf.h>
#include <patch/affine.h>
#include <vote.h>
#include <texture.h>
#include <map>

template <typename Patch, typename Scalar>
void typedMexFunction(int nout, mxArray *out[], int nin, const mxArray *in[],
		const pm::Texture &source, const pm::Texture &target, pm::VoteParams &params);

#define mexArgs nout, out, nin, in

/**
 * Usage:
 * 
 * [vote, weights] = votemex( source, target, nnf, options )
 * 
 */
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
	
	// checking the input
	if (nin != 4) {
		mexErrMsgIdAndTxt("MATLAB:vote:invalidNumInputs",
				"Required arguments: source, target, nnf1, [nnf2,] options.");
	}
	// checking the output
	if (nout > 2) {
		mexErrMsgIdAndTxt("MATLAB:vote:maxlhs",
				"Too many output arguments.");
	}
	// checking the arguments
	// - loading the images
	pm::Texture source, target;
	mxLoadTexture(source, in[0], "Corrupted source!");
	mxLoadTexture(target, in[1], "Corrupted target!");

	// options for the patch size
	const mxArray *options = in[3];
	if (!mxIsStruct(options)) {
		mexErrMsgIdAndTxt("MATLAB:vote:options",
				"The option parameter was not found!");
	}
	mxArray *tmp;
	if (tmp = mxGetField(options, 0, "num_threads")) {
		int num = mxCheckedScalar(tmp, "Invalid num_threads");
		 omp_set_num_threads(num);
	}
	if (tmp = mxGetField(options, 0, "rand_seed")) {
		unsigned int s = mxCheckedScalar(tmp, "Invalid seed number!");
		pm::seed(s);
	} else {
		pm::seed();
	}
	if (tmp = mxGetField(options, 0, "interpolation")) {
		if(mxStringEquals(tmp, "nn"))
			source.interp = target.interp = pm::NearestNeighbor;
		else if(mxStringEquals(tmp, "bilinear"))
			source.interp = target.interp = pm::Bilinear;
	}
	if (tmp = mxGetField(options, 0, "boundaries")) {
		if(mxStringEquals(tmp, "unchecked"))
			source.model = target.model = pm::Unchecked;
		else if(mxStringEquals(tmp, "checked"))
			source.model = target.model = pm::Checked;
		else if(mxStringEquals(tmp, "same"))
			source.model = target.model = pm::Same;
		else if(mxStringEquals(tmp, "toroidal"))
			source.model = target.model = pm::Toroidal;
	}
    
    // channel mapping
    std::map<int, int> channel_map;
    if (tmp = mxGetField(options, 0, "vote_channels")) {
        std::vector<int> vote_channels;
        mxLoadVector(vote_channels, tmp, "Invalid vote_channels!");
        for(int c_i = 0; c_i < vote_channels.size(); ++c_i) {
            channel_map[vote_channels[c_i]] = c_i;
        }
    }

	// vote parameters
	pm::VoteParams params;
	if (tmp = mxGetField(options, 0, "patch_size")) {
		params.patchSize = mxCheckedScalar(tmp, "Invalid patch size number.");
	}
	if (tmp = mxGetField(options, 0, "vote_method")) {
		if (!mxIsChar(tmp)) {
			mexErrMsgIdAndTxt("MATLAB:vote:wrongMethod",
					"Invalid voting method!");
		}
		// DEFAULT /////////////////////////////////////////////////////////////
		if (mxStringEquals(tmp, "default")) 
			params.method = pm::DefaultVoting;
		// BIDIR ///////////////////////////////////////////////////////////////
		else if (mxStringEquals(tmp, "bidir_sim")){
			params.method = pm::BiDirSimVoting;
			// relative weight of reverse nnf
			if(tmp = mxGetField(options, 0, "reverse_nnf_weight"))
				params.bidirSimWeight = mxCheckedScalar(tmp, "Invalid reverse_nnf_weight.");
		} 
		// GBA /////////////////////////////////////////////////////////////////
		else if (mxStringEquals(tmp, "default_gba")) {
			params.method = pm::DefaultGBAVoting;
			mxLoadGB<pm::GainBias<float> >(options);
		}
		// DIRECT //////////////////////////////////////////////////////////////
		else if(mxStringEquals(tmp, "direct"))
			params.method = pm::DirectVoting;
		// MEDIAN //////////////////////////////////////////////////////////////
		else if(mxStringEquals(tmp, "median"))
			params.method = pm::MedianVoting;
		// TILED ///////////////////////////////////////////////////////////////
		else if (mxStringEquals(tmp, "tiled")) 
			params.method = pm::TiledVoting;
		// MEAN-SHIFT //////////////////////////////////////////////////////////
		else if (mxStringEquals(tmp, "meanshift")) {
			params.method = pm::MeanShiftVoting;
			params.meanShift.windows.resize(source.channels());
			if (tmp = mxGetField(options, 0, "meanshift_window")) {
				mxLoadScalars(&params.meanShift.windows[0], source.channels(), tmp, "Invalid mean-shift window!");
			}
			if (tmp = mxGetField(options, 0, "meanshift_threshold")) {
				if(mxIsNumeric(tmp))
					params.meanShift.threshold = mxGetScalar(tmp);
				else
					mexErrMsgIdAndTxt("MATLAB:vote:meanshiftThreshold",
							"The mean-shift threshold must be a scalar!");
			}
			if (tmp = mxGetField(options, 0, "meanshift_shuffle")) {
				params.meanShift.shuffle = mxCheckedScalar(tmp, "meanshift_shuffle must be a number!") != 0;
			}
		} 
		// HISTOGRAM ///////////////////////////////////////////////////////////
		else if(mxStringEquals(tmp, "histogram")) { 
			params.method = pm::HistogramVoting;
			if(mxHasField(options, 0, "reverse_nnf")) {
				// we use bidirectional similarity with histogram voting
				params.method = pm::BiDirSimWithHistogramVoting;
				// relative weight of reverse nnf
				if(tmp = mxGetField(options, 0, "reverse_nnf_weight"))
					params.bidirSimWeight = mxCheckedScalar(tmp, "Invalid reverse_nnf_weight.");
			} else if(tmp = mxGetField(options, 0, "feature_method")) {
				params.method = pm::FeatureWithHistogramVoting;
				if(mxStringEquals(tmp, "min")){
					params.featureOp = &std::min<float>;
					params.featureDefault = 1e10;
				}else if(mxStringEquals(tmp, "max")){
					params.featureOp = &std::max<float>;
					params.featureDefault = -1e10;
				} else mexErrMsgIdAndTxt("MATLAB:vote:featureOp", "Invalid feature_method");
			}
            params.histNormalize = mxBoolField(options, 0, "hist_normalize", true);
            if (tmp = mxGetField(options, 0, "hist_params")) {
                // format
                //      channel bins min max weight [boost]
                std::vector< double > hp;
                mxLoadVector(hp, tmp, "Invalid hist_params");
                const int h = mxGetM(tmp);
                const int w = mxGetN(tmp);
                for(int p = 0; p < h; ++p) {
                    params.histChannels.push_back(hp[p + 0 * h]);
                    params.histBins.push_back(hp[p + 1 * h]);
                    params.histRanges.push_back(pm::Range(hp[p + 2 * h], hp[p + 3 * h]));
                    params.histWeights.push_back(hp[p + 4 * h]);
                    if(w == 6) params.histBoosts.push_back(hp[p + 5 * h]);
                }
            }
            // check the sizes
            int N = params.histBins.size();
            if(N != params.histChannels.size() || N != params.histRanges.size() || N != params.histWeights.size()){
                mexErrMsgIdAndTxt("MATLAB:vote:hist_num", "The histogram parameters don't have the same size!");
            }
            // convert the channels
            for (int c = 0; c < params.histChannels.size(); ++c) {
                int c_i = params.histChannels[c];
                std::map<int, int>::const_iterator it = channel_map.find(c_i);
                if(it != channel_map.end()) {
                    c_i = it->second; // convert using the map from vote_channels
                }
				// std::cout << "Hist: mapping " << params.histChannels[c] << " to " << c_i << "\n";
                if(c_i < 0 || c_i >= source.channels()) {
                    mexErrMsgIdAndTxt("MATLAB:vote:hist_channels", "Trying to access channel #%d of %d!", c_i + 1, source.channels());
                }
				// reset back
				params.histChannels[c] = c_i;
            }
			// default uses a histogram for the L channel
			if(N == 0) {
				params.histChannels.push_back(0);
				params.histBins.push_back(100);
				params.histRanges.push_back(pm::Range(0.0f, 100.0f));
				params.histWeights.push_back(0.1f);
			}
			for(int c = 0; c < params.histChannels.size(); ++c) {
				int c_i = params.histChannels[c];
				if(c_i < 0 || c_i >= source.channels()) {
                    mexErrMsgIdAndTxt("MATLAB:vote:hist_channels", "Trying to access channel #%d of %d!", c_i + 1, source.channels());
                }
			}
		}
		// WEIGHTED ////////////////////////////////////////////////////////////
		else if (mxStringEquals(tmp, "weighted")) {
			params.method = pm::WeightedVoting;
			if((tmp = mxGetField(options, 0, "vote_mask")) && mxGetNumberOfElements(tmp) > 1) {
				params.weightMask = mxArrayToImage(tmp, "Invalid weight mask");
			} else
				mexErrMsgIdAndTxt("MATLAB:vote:vote_mask", "Missing the vote weight mask.");
            if (tmp = mxGetField(options, 0, "vote_base"))
                params.weightBase = mxCheckedScalar(tmp, "Invalid vote_base parameter.");
		}
		// INVALID /////////////////////////////////////////////////////////////
		else {
			mexErrMsgIdAndTxt("MATLAB:vote:wrongMethod",
					"Valid voting methods include: "
					"default, tiled, gba, meanshift, histogram");
		}
	}
	params.filter = pm::Filter(params.patchSize);
    if (tmp = mxGetField(options, 0, "binary_channels")) {
        mxLoadVector(params.binaryChannels, tmp, "Invalid hist_params");
        for(int i = 0; i < params.binaryChannels.size(); ++i) {
            int c_i = params.binaryChannels[i];
            std::map<int, int>::const_iterator it = channel_map.find(c_i);
            if(it != channel_map.end()) {
                c_i = it->second; // convert using the map from vote_channels
            } else {
                c_i = i;
            }
            if(c_i < 0 || c_i >= source.channels()) {
                mexErrMsgIdAndTxt("MATLAB:vote:binary_channels", "Trying to access channel #%d of %d!", c_i + 1, source.channels());
            }
            params.binaryChannels[i] = c_i;
        }
    }
	if (tmp = mxGetField(options, 0, "vote_weight")) {
		mxLoadFilter(params.filter, tmp, "Invalid vote_weight filter!");
	}
	// do we want the weights out?
	if (nout >= 2) {
		if(tmp = mxGetField(options, 0, "weight_type")){
			if(mxStringEquals(tmp, "variance"))
				params.weightType = pm::PixelVariance;
			else
				params.weightType = pm::PixelWeight;
		} else {
			params.weightType = pm::PixelWeight;
		}
	}

	pm::PatchType patchType = pm::BasicFloat;
	if (tmp = mxGetField(options, 0, "patch_type")){
		if (mxStringEquals(tmp, "affine"))
			patchType = pm::Affine;
		else if(mxStringEquals(tmp, "basicf"))
			patchType = pm::BasicFloat;
		else if(mxStringEquals(tmp, "basic"))
			patchType = pm::Basic;
	}
	
	switch(patchType) {
		case pm::Basic: typedMexFunction<pm::BasicPatchX, float>(
				mexArgs, source, target, params);
			break;
		case pm::BasicFloat: typedMexFunction<pm::BasicFloatPatchX, float>(
				mexArgs, source, target, params);
			break;
		case pm::Affine: typedMexFunction<pm::AffinePatchX, float>(
				mexArgs, source, target, params);
			break;
		default:
			mexErrMsgIdAndTxt("MATLAB:nnf:patch_type", "Invalid patch type!");
	}
	
}

template <typename Patch, typename Scalar>
void typedMexFunction(int nout, mxArray *out[], int nin, const mxArray *in[],
		const pm::Texture &source, const pm::Texture &target, pm::VoteParams &params) {
	
	// - loading the nnf map
	Patch::width(params.patchSize);
	typedef pm::NearestNeighborField<Patch, Scalar> NNF;
	
	if(params.weightType != pm::NoWeight) {
		params.weights = new float[source.rows * source.cols]; // no need to initialize
	}
	
	// voting
	Image result;
	if (params.method == pm::BiDirSimVoting || params.method == pm::BiDirSimWithHistogramVoting){
		// using the T to S and the S to T nnfs
		// /!\ uses placement-new initialization for our array since def constructor
		//	   is not the good choice here
		// @see http://stackoverflow.com/questions/4754763/c-object-array-initialization-without-default-constructor
		NNF *nnf = static_cast<NNF *>(operator new[](2 * sizeof(NNF)));
		new(&nnf[0])NNF(&source, &target, false);
		new(&nnf[1])NNF(&target, &source, false);
		// load both nnfs' data
		mxLoadNNF(nnf[0], in[2], "Third argument should be a valid nnf map!");
		if(const mxArray *revNNF = mxGetField(in[3], 0, "reverse_nnf")){
			mxLoadNNF(nnf[1], revNNF, "Bidirectional voting requires a valid reverse_nnf option!");
		} else {
			mexErrMsgIdAndTxt("MATLAB:vote:reverse_nnf", "Missing a reverse_nnf for bidirectional voting!");
		}
		
		// vote with both nnfs
		result = pm::vote(&nnf[0], params);
		
		// delete the nnf data (reverse order)
		// destruct in inverse order
		nnf[1].~NNF();
		nnf[0].~NNF();
		operator delete[]( static_cast<void *>(nnf) );
	} else {
		NNF nnf(&source, &target, false); // no need for the distances
		mxLoadNNF(nnf, in[2], "Third argument should be a valid nnf map!");
		// using the T to S nnf
		result = pm::vote(&nnf, params);
	}

	// output
	out[0] = mxImageToArray(result);
	if (nout >= 2) {
		out[1] = mxCreateMatrix(result.rows, result.cols, mxSINGLE_CLASS);
		Scalar *data = reinterpret_cast<Scalar *> (mxGetData(out[1]));
		float *weights = params.weights;
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
		for (int y = 0; y < result.rows; ++y) {
			for (int x = 0; x < result.cols; ++x) {
				data[y + x * result.rows] = weights[result.cols * y + x]; // transposed!
			}
		}
		// free it since we took ownership!
		delete[] weights;
	}
}
