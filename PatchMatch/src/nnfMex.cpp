/*******************************************************************************
 * nnfMex.cpp - nnf matlab function implementation
 *******************************************************************************
 * Add license here...
 *******************************/

// other stuff
#include "mexwrap.h"
// our stuff
#include <distance.h>
#include <gb.h>
#include <gpm.h>
#include <patch/basic.h>
#include <patch/affine.h>
#include <patch/multires.h>
#include <texture.h>

#define mexArgs nout, out, nin, in
#define SOURCE_IDX	0
#define TARGET_IDX	1
#define NNF_IDX		2
#define MASK_IDX	3
#define OPTIONS_IDX	4

template <typename Patch, typename Scalar>
void typedMexFunction(int nout, mxArray *out[], int nin, const mxArray *in[],
		const pm::Texture *source, const pm::Texture *target, pm::NNSettings &settings);

/**
 * Usage:
 * 
 * [newNNF, dist, conv] = nnfmex( source, target, prevNNF, mask, options )
 * 
 */
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
	
	// checking the input
	if (nin != 5) {
		mexErrMsgIdAndTxt("MATLAB:nnf:invalidNumInputs",
				"Requires 5 arguments! (#in = %d)", nin);
	}
	// checking the output
	if (nout > 3) {
		mexErrMsgIdAndTxt("MATLAB:nnf:maxlhs",
				"Too many output arguments.");
	}
	
	// options parameter
	if (!mxIsStruct(in[OPTIONS_IDX])) {
		mexErrMsgIdAndTxt("MATLAB:nnf:invalidOptions",
				"Options must be put in a structure.");
	}
	const mxArray *options = in[OPTIONS_IDX];
	
	// checking the arguments
	pm::Texture source, *target;
	mxLoadTexture(source, in[SOURCE_IDX], "Corrupted source!");
	pm::Texture localTarget;
	if(mxBoolField(options, 0, "self_nnf")) {
		target = &source;
	} else {
		mxLoadTexture(localTarget, in[TARGET_IDX], "Corrupted target!");
		target = &localTarget;
	}

	// arguments
	pm::NNSettings settings;
	settings.storeConvergence = nout >= 3;
	mxArray *tmp;
	// debug purpose
	if (tmp = mxGetField(options, 0, "segfault")) {
		float *ptr = NULL; *ptr = 1;
	}
	if (tmp = mxGetField(options, 0, "num_threads")) {
		int num = mxCheckedScalar(tmp, "Invalid num_threads");
		 omp_set_num_threads(num);
	}
	// --- iterations ----------------------------------------------------------
	if (tmp = mxGetField(options, 0, "iterations")) {
		settings.iterations = mxCheckedScalar(tmp, "Invalid iteration number.");
	}
	if (tmp = mxGetField(options, 0, "it_until_conv")) {
		settings.untilConvergence = mxCheckedScalar(tmp, "Invalid it_until_conv parameter!");
	}
	if (tmp = mxGetField(options, 0, "it_start_order")) {
		settings.startOrder = mxCheckedScalar(tmp, "Invalid it_start_order!") == 0 ? pm::Normal : pm::Reverse;
	}
	// --- random search -------------------------------------------------------
	if (tmp = mxGetField(options, 0, "coh_rand_search")) {
		settings.coherentRandSearch = mxCheckedScalar(tmp, "Invalid coherent random search parameter.") != 0;
	}
	if (tmp = mxGetField(options, 0, "rand_search")) {
		settings.randSearch = mxCheckedScalar(tmp, "Invalid rand_search!");
		if (tmp = mxGetField(options, 0, "max_rand_search")){
			settings.maxRandSearch = mxCheckedScalar(tmp, "Invalid max_rand_search!");
		}
	}
	if (tmp = mxGetField(options, 0, "mixed_rand_search")){
		settings.mixedRandSearch = mxCheckedScalar(tmp, "Invalid mixed_rand_search!") != 0;
	}
	if (tmp = mxGetField(options, 0, "window_size")) {
		settings.windowSize = mxCheckedScalar(tmp, "Invalid window size number.");
	}
	if (tmp = mxGetField(options, 0, "window_size_decrease")) {
		settings.windowSizeDecr = mxCheckedScalar(tmp, "Invalid window_size_decrease parameter.") != 0;
	}
	// --- aligned search ------------------------------------------------------
	if(tmp = mxGetField(options, 0, "aligned_search")) {
		settings.alignedSearch = mxCheckedScalar(tmp, "Invalid aligned_search parameter.");
		if(settings.alignedSearch > 0) {
			if(tmp = mxGetField(options, 0, "generators")) {
				pm::MatXD G(tmp);
				settings.alignG1 = pm::Point2f(G.read<float>(0, 1), G.read<float>(0, 0));
				settings.alignG2 = pm::Point2f(G.read<float>(1, 1), G.read<float>(1, 0));
				settings.alignP1 = G.read<float>(0, 2);
				settings.alignP2 = G.read<float>(1, 2);
			} else {
				mexErrMsgIdAndTxt("MATLAB:nnf:aligned_search", "Had aligned_search, but no generator!");
			}
			if(tmp = mxGetField(options, 0, "align_jitter")){
				settings.alignJitter = mxCheckedScalar(tmp, "Invalid align_jitter parameter!");
			}
			if (tmp = mxGetField(options, 0, "max_aligned_search")){
				settings.maxAlignedSearch = mxCheckedScalar(tmp, "Invalid max_aligned_search!");
			}
		}
	}
	// --- nnf mask ------------------------------------------------------------
	if (mxGetNumberOfElements(in[MASK_IDX]) > 1) {
		settings.mask = mxArrayToImage(in[MASK_IDX], "Invalid mask");
		if (tmp = mxGetField(options, 0, "scramble_nnf")) {
			settings.scramble = mxCheckedScalar(tmp, "Invalid scramble_nnf parameter.") != 0;
		}
	}
	// --- completeness --------------------------------------------------------
	if ((settings.completeness.penalty = mxScalarField(options, 0, "comp_penalty")) > 0) {
		if(tmp = mxGetField(options, 0, "comp_weight")){
			settings.completeness.weight = mxCheckedScalar(tmp, "Invalid completeness weight.");
			settings.completeness.invWeight = 1.0f / settings.completeness.weight;
		}
		if(tmp = mxGetField(options, 0, "comp_exponent")) 
			settings.completeness.exponent = mxCheckedScalar(tmp, "Invalid completeness exponent.");
		if(tmp = mxGetField(options, 0, "comp_abs_threshold")) 
			settings.completeness.absoluteThreshold = mxCheckedScalar(tmp, "Invalid completeness absolute threshold");
		if(tmp = mxGetField(options, 0, "comp_rel_threshold")) 
			settings.completeness.relativeThreshold = mxCheckedScalar(tmp, "Invalid completeness relative threshold");
        if(tmp = mxGetField(options, 0, "comp_dist"))
            settings.occDist = mxCheckedScalar(tmp, "Invalid completeness distance boolean.") != 0;
        settings.storeOccRatio = mxBoolField(options, 0, "saveOccRatio", false);
		settings.completePropagation = mxBoolField(options, 0, "comp_prop", false);
		
		// --- incomplete search -----------------------------------------------
		if((settings.incompleteSearch = mxScalarField(options, 0, "incomp_search")) > 0) {
			if(tmp = mxGetField(options, 0, "max_incomp_search"))
				settings.maxIncompleteSearch = mxCheckedScalar(tmp, "Invalid max_incomp_search");
			if(tmp = mxGetField(options, 0, "incomp_sigma_ratio"))
				settings.jumpSigmaRatio = mxCheckedScalar(tmp, "Invalid incomp_sigma_ratio");
			if(tmp = mxGetField(options, 0, "incomp_buffer_size"))
				settings.jumpBufferSize = mxCheckedScalar(tmp, "Invalid incomp_buffer_size");
		}
	}
	// --- patch type ----------------------------------------------------------
	if (tmp = mxGetField(options, 0, "patch_size")) {
		settings.patchSize = mxCheckedScalar(tmp, "Invalid patch size number.");
	}
	if (tmp = mxGetField(options, 0, "multires")) {
		settings.multires = mxCheckedScalar(tmp, "Invalid multires depth.");
	}
	// --- texture model -------------------------------------------------------
	if (tmp = mxGetField(options, 0, "interpolation")) {
		if(mxStringEquals(tmp, "nn"))
			source.interp = localTarget.interp = pm::NearestNeighbor;
		else if(mxStringEquals(tmp, "bilinear"))
			source.interp = localTarget.interp = pm::Bilinear;
	}
	if (tmp = mxGetField(options, 0, "boundaries")) {
		if(mxStringEquals(tmp, "unchecked"))
			source.model = localTarget.model = pm::Unchecked;
		else if(mxStringEquals(tmp, "checked"))
			source.model = localTarget.model = pm::Checked;
		else if(mxStringEquals(tmp, "same"))
			source.model = localTarget.model = pm::Same;
		else if(mxStringEquals(tmp, "toroidal"))
			source.model = localTarget.model = pm::Toroidal;
	}
	// --- random seed ---------------------------------------------------------
	if (tmp = mxGetField(options, 0, "rand_seed")) {
		unsigned int s = mxCheckedScalar(tmp, "Invalid seed number!");
		pm::seed(s);
	} else {
		pm::seed();
	}
	// --- distance type -------------------------------------------------------
	if (tmp = mxGetField(options, 0, "dist_type")) {
		if (mxStringEquals(tmp, "gbassd")){
			settings.distType = pm::GBASSD;
			mxLoadGB<pm::GainBias<float> >(options); // g+b parameters
		} else if (mxStringEquals(tmp, "ssd")){
			settings.distType = pm::SSD;
        } else if (mxStringEquals(tmp, "lp")){
            settings.distType = pm::LP;
            if(tmp = mxGetField(options, 0, "p_space")){
                float p = mxCheckedScalar(tmp, "Invalid p_space value");
                if(p <= 0.0f)
                    mexErrMsgIdAndTxt("MATLAB:nnf:p_space", "p_space must be strictly positive!");
                pm::dist::pSpace(p);
            }
            pm::dist::pSpaceRoot(mxBoolField(options, 0, "p_root", true));
		} else if (mxStringEquals(tmp, "wssd")) {
			settings.distType = pm::wSSD;
			if(tmp = mxGetField(options, 0, "dist_weight")) {
				pm::Texture &w = pm::dist::ssdWeight();
				mxLoadTexture(w, tmp, "Invalid distance weight");
				w.model = pm::Same;
				// we check it
				int nnz = 0;
				for(int y = 0; y < w.height; ++y) {
					for(int x = 0; x < w.width; ++x) {
						const pm::Texture &weights = pm::dist::ssdWeight();
						if(weights.at<float>(y, x) > 0.0f) {
							++nnz;
						}
					}
				}
				std::cout << "wssd with " << nnz << " non-zero weights";
			}
		}
		//else if (mxStringEquals(tmp, "parssd"))
		//	settings.distType = pm::ParSSD;
	}
	// --- identity distance
	if (tmp = mxGetField(options, 0, "min_patch_disp")){
		settings.minPatchDisp = mxCheckedScalar(tmp, "Invalid min_patch_disp.");
	}
    // added by huajie. weight sequence dataset, represent the index of weight
    if (tmp = mxGetField(options, 0, "weight")){
		settings.weightIndex = mxCheckedScalar(tmp, "Invalid weightIndex.");
	}
	// type dependent nnf
	pm::PatchType patchType = settings.multires <= 0 ? pm::BasicFloat : pm::BasicFloatMultiRes;
	if (tmp = mxGetField(options, 0, "patch_type")){
		if (mxStringEquals(tmp, "basic"))
			patchType = settings.multires <= 0 ? pm::Basic : pm::BasicMultiRes;
		else if (mxStringEquals(tmp, "basicf"))
			patchType = settings.multires <= 0 ? pm::BasicFloat : pm::BasicFloatMultiRes;
		else if (mxStringEquals(tmp, "affine"))
			patchType = settings.multires <= 0 ? pm::Affine : pm::AffineMultiRes;
		// else if(mxStringEquals(tmp, "bilinear_affine"))
		// 	patchType = pm::BilinearAffine;
		else 
			mexErrMsgIdAndTxt("MATLAB:nnf:patch_type", "Invalid patch type!");
	}
	
#define CALL(patch, scalar) typedMexFunction<patch, scalar>(mexArgs, &source, target, settings)
	
	switch(patchType) {
		case pm::Basic: CALL(pm::BasicPatchX, float); break;
		case pm::BasicFloat: CALL(pm::BasicFloatPatchX, float); break;
		case pm::Affine: CALL(pm::AffinePatchX, float); break;
		case pm::BasicMultiRes: CALL(pm::MultiRes<pm::BasicPatchX>, float); break;
		case pm::BasicFloatMultiRes: CALL(pm::MultiRes<pm::BasicFloatPatchX>, float); break;
		case pm::AffineMultiRes: CALL(pm::MultiRes<pm::AffinePatchX>, float); break;
		default:
			mexErrMsgIdAndTxt("MATLAB:nnf:patch_type", "Invalid patch type!");
	}
}


template <typename Patch, typename Scalar>
void typedMexFunction(int nout, mxArray *out[], int nin, const mxArray *in[],
		const pm::Texture *source, const pm::Texture *target, pm::NNSettings &settings){
	
	Patch::width(settings.patchSize); // dynamic size
	if(settings.multires > 0){
		Patch::depth(settings.multires); // multires depth
		std::cout << "Multires: " << Patch::depth() << "\n";
	}
	typedef pm::NearestNeighborField<Patch, Scalar> NNF;

	// previous NNF
	NNF *prevNNF = NULL;
	if (mxGetNumberOfElements(in[NNF_IDX]) > 1) {
		prevNNF = new NNF(source, target);
		mxLoadNNF(*prevNNF, in[NNF_IDX], "Invalid prevNNF parameter!");
	}
	
	const mxArray *tmp;
	typedef typename NNF::Extension ExtNNF;
	ExtNNF *extNNF = NULL;
	if(tmp = mxGetField(in[OPTIONS_IDX], 0, "ext_nnf")){
		extNNF = new ExtNNF(target, target);
		mxLoadNNF(*extNNF, tmp, "Invalid extNNF parameter");
		
	}

	const NNF *nnf = pm::nnf<Patch, Scalar>(source, target, prevNNF, extNNF, settings);
	if (nnf == NULL) {
		mexErrMsgIdAndTxt("MATLAB:nnf:zeronnf",
				"NNF is null!");
	}

	// output
	// 1 = nnf
	const int nnfDims = Patch::dimensions();
	Image nnfImg(nnf->height, nnf->width, IM_MAKETYPE(DataDepth<Scalar>::value, nnfDims));
	for (int y = 0; y < nnf->height; ++y) {
		for (int x = 0; x < nnf->width; ++x) {
			Scalar *dptr = nnfImg.ptr<Scalar>(y, x);
			const Patch &p = nnf->get(y, x);
			p.store(dptr, nnfDims);
		}
	}
	out[0] = mxImageToArray(nnfImg);

	// 2 = distance
	if (nout >= 2) {
		pm::MatXD dist(nnf->height, nnf->width, classOf(*source));
		for (int y = 0; y < nnf->height; ++y) {
			for (int x = 0; x < nnf->width; ++x) {
				dist.update(y, x, settings.occDist ? nnf->occDensity(y, x) : nnf->distance(y, x));
			}
		}
		out[1] = dist;
	}
	
	// 3 = convergence data
	if (nout == 3) {
		const char *fields[] = {
			"prop_count", "rs_count", "as_count", "is_count", 
			"max_dist", "mean_dist", "coh_ratio", "occ_ratio"
		};
        const unsigned int F = settings.storeOccRatio ? 8 : 7;
		out[2] = mxCreateStructMatrix(1, 1, F, fields);
		int N = settings.convergence.size();
		std::cout << "N=" << N << "\n";
		float *ptr[8];
		for(int f = 0; f < F; ++f) {
			mxArray *arr = mxCreateMatrix(1, N);
			mxSetFieldByNumber(out[2], 0, f, arr);
			ptr[f] = reinterpret_cast<float *>(mxGetData(arr));
		}
		// fill these arrays
		for(int i = 0; i < settings.convergence.size(); ++i) {
			pm::ConvergenceData &cd = settings.convergence[i];
			ptr[0][i] = cd.propCount;
			ptr[1][i] = cd.rsCount;
            ptr[2][i] = cd.asCount;
			ptr[3][i] = cd.isCount;
			ptr[4][i] = cd.maxDist;
			ptr[5][i] = cd.meanDist;
			ptr[6][i] = cd.cohRatio;
            if(settings.storeOccRatio)
                ptr[7][i] = cd.occRatio;
		}
	}
	// free all the nnf data
	if (nnf != NULL) delete nnf;
	// Note: no need to delete the previous nnf
	//	as it will be the output instance (reused by pm::nnf)
	// => do not do: if(prevNNF != NULL) delete prevNNF;
}