/*******************************************************************************
 * segmentMex.cpp - nnf segmentation for matlab
 *******************************************************************************
 * Add license here...
 *******************************/

// other stuff
#include <mexutil.h>
#include <boost/unordered_map.hpp>
// our stuff
#include <nnf.h>
#include <patch/affine.h>
#include <texture.h>
#include <segment.h>

template <typename Patch, typename Scalar>
void typedMexFunction(int nout, mxArray *out[], const mxArray *nnfArray);

#define mexArgs nout, out, in[0]

/**
 * Usage:
 * 
 * [boundaries, labelmap] = segmentmex( nnf, patch_type )
 * 
 */
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
	
	// checking the input
	if (nin > 2) {
		mexErrMsgIdAndTxt("MATLAB:segment:invalidNumInputs",
				"Only two arguments allowed.");
	}
	// checking the output
	if (nout > 2) {
		mexErrMsgIdAndTxt("MATLAB:segment:maxlhs",
				"Too many output arguments.");
	}

	pm::PatchType patchType = pm::Basic;
	if (nin == 2){
        const mxArray *patch_type = in[1];
		if (mxStringEquals(patch_type, "affine"))
			patchType = pm::Affine;
		else if(mxStringEquals(patch_type, "basicf"))
			patchType = pm::BasicFloat;
		else if(mxStringEquals(patch_type, "basic"))
			patchType = pm::Basic;
	}
	
	switch(patchType) {
		case pm::Basic: typedMexFunction<pm::BasicPatchX, float>(mexArgs); break;
		case pm::BasicFloat: typedMexFunction<pm::BasicFloatPatchX, float>(mexArgs); break;
		case pm::Affine: typedMexFunction<pm::AffinePatchX, float>(mexArgs); break;
		default:
			mexErrMsgIdAndTxt("MATLAB:nnf:patch_type", "Invalid patch type!");
	}
	
}

template <typename Patch, typename Scalar>
void typedMexFunction(int nout, mxArray *out[], const mxArray *nnfArray) {
	
    // - size
    const int h = mxGetDimensions(nnfArray)[0];
	const int w = mxGetDimensions(nnfArray)[1];
    
	// - loading the nnf map
	Patch::width(10);
	typedef pm::NearestNeighborField<Patch, Scalar> NNF;
	NNF nnf(w, h); // no need for the distances
	mxLoadNNF(nnf, nnfArray, "Invalid nnf!");
    
    // - segmentation
    pm::BaseSegmentation* s = pm::segment<pm::NoData>(&nnf);

	// - boundaries
	switch(nout){
		case 2:
		{
			// label storage
			out[1] = mxCreateMatrix(h, w, mxSINGLE_CLASS);
			float *data = reinterpret_cast<float *>(mxGetData(out[1])); // unsigned int?
			// label map
			// create simplest labels from 1 to N
			unsigned int N = 1;
			typedef pm::BaseSegmentation::Label Label;
			typedef boost::unordered_map<Label, unsigned int> LabelMap;
			LabelMap labelMap;
			// loop
			for (int y = 0; y < s->rows; ++y) {
				for (int x = 0; x < s->cols; ++x) {
					// map the label to a simple one
					unsigned int simpleLabel; //< to find
					Label label = s->at(y, x).label();
					LabelMap::iterator it = labelMap.find(label);
					if(it == labelMap.end()) {
						simpleLabel = N++; // new label
						labelMap[label] = simpleLabel;
					} else
						simpleLabel = it->second; // already used before
					
					// store the simple label
					data[y + x * s->rows] = simpleLabel;
				}
			}
		}
		case 1:
		default:
		{
			out[0] = mxCreateMatrix(h, w, mxLOGICAL_CLASS);
			bool *data = reinterpret_cast<bool *>(mxGetData(out[0]));
#if _OPENMP
#pragma omp parallel for collapse(2)
#endif
			for (int y = 0; y < s->rows; ++y) {
				for (int x = 0; x < s->cols; ++x) {
					data[y + x * s->rows] = s->at(y, x).isBoundary();
				}
			}
		}
	}
	// release memory
	if (s != NULL) delete s;
}

