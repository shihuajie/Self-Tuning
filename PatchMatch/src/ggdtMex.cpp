/*******************************************************************************
 * ggdtMex.cpp - generalized geodesic distance transform for matlab
 *******************************************************************************
 * Add license here...
 *******************************/

#include "ggdt.h"
#include "mexwrap.h"

typedef pm::MatWrapper FastImage;

/**
 * Usage:
 * 
 * D = ggdtmex( M, mu, I, delta, num_iters )
 * 
 */
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
	// checking the input
	if (nin > 5) {
		mexErrMsgIdAndTxt("MATLAB:ggdt:maxrhs",
				"Too many arguments, max 5, got %d", nin);
	} else if(nin < 2) {
		mexErrMsgIdAndTxt("MATLAB:ggdt:minrhs",
				"Requires at least two arguments: M and mu");
	}
	
	// checking the output
	if (nout > 1) {
		mexErrMsgIdAndTxt("MATLAB:ggdt:maxlhs", "Too many output arguments.");
	} else if(nout < 1) {
		mexErrMsgIdAndTxt("MATLAB:ggdt:minlhs", "Requires an output argument.");
	}
	
	// read the input gray mask
	const FastImage M(in[0]);
	FastImage I;
	double mu = mxCheckedScalar(in[1], "Invalid mu!");
	double delta = 0.0;
	int iterations = 5;
	switch(nin) {
		case 5:
			I = FastImage(in[mxGetNumberOfElements(in[2]) > 1 ? 2 : 0]);
			delta = mxCheckedScalar(in[3], "Invalid delta!");
			iterations = mxCheckedScalar(in[4], "Invalid number of iterations!");
			break;
		case 4:
			if(mxGetNumberOfElements(in[2]) > 1) {
				I = FastImage(in[2]);
				delta = 1.0;
			} else {
				I = FastImage(in[0]);
				delta = mxCheckedScalar(in[2], "Invalid delta!");
			}
			iterations = mxCheckedScalar(in[3], "Invalid number of iterations!");
		case 3:
			delta = 0.0;
			break;
	}
	
	switch(M.depth()) {
		case IM_32F:
			out[0] = pm::ggdt<float>(M, mu, &I, float(delta), iterations);
			break;
		case IM_64F:
			out[0] = pm::ggdt<double>(M, mu, &I, delta, iterations);
			break;
		default:
			mexErrMsgIdAndTxt("MATLAB:ggdt:class", "Only single or double supported!");
	}
}