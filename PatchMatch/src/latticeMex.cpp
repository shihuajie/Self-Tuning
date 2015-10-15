/*******************************************************************************
 * latticeMex.cpp - lattice extraction for matlab
 *******************************************************************************
 * Add license here...
 *******************************/

#include "lattice.h"
#include "mexwrap.h"

typedef pm::MatWrapper MatXD;

/**
 * Usage:
 * 
 * [ [G1;G2], [p1;p2] ] = ggdtmex( X, options )
 * [ map, [minY,minX] ]	= ggdtmex( X, 'map', options )
 * [ E ]				= ggdtmex( X, 'energy', options )
 * 
 */
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[]) {
	// checking the input
	if (nin > 3) {
		mexErrMsgIdAndTxt("MATLAB:lattice:maxrhs",
				"Too many arguments, max 3, got %d", nin);
	} else if (nin < 1) {
		mexErrMsgIdAndTxt("MATLAB:lattice:minrhs",
				"Requires at least the offset matrix");
	}

	// checking the output
	if (nout > 2) {
		mexErrMsgIdAndTxt("MATLAB:lattice:maxlhs", "Too many output arguments.");
	} else if (nout < 1) {
		mexErrMsgIdAndTxt("MATLAB:lattice:minlhs", "Requires an output argument.");
	}

	// testing the input offsets X (Nx2)
	const MatXD X(in[0]);
	if (X.rows == 0) mexErrMsgIdAndTxt("MATLAB:lattice:X", "Empty offsets!");
	if (X.cols != 2) mexErrMsgIdAndTxt("MATLAB:lattice:X", "Invalid offset size, should be Nx2!");

	// loading X
	std::vector<pm::Point2f> points;
	points.reserve(X.rows);
	for (unsigned int i = 0; i < X.rows; ++i) {
		points.push_back(pm::Point2f(
				X.read<float>(i, 0),
				X.read<float>(i, 1))
				);
	}
	
	// mode
	enum {
		Energy,
		Map,
		Result
	} mode = Result;
	const mxArray *options = NULL;
	if (nin >= 2) {
		if (mxIsChar(in[1])) {
			if(nin >= 3) {
				options = in[2];
			}
			if(mxStringEquals(in[1], "map"))
				mode = Map;
			else if(mxStringEquals(in[1], "energy")){
				mode = Energy;
			}
		} else {
			options = in[1];
		}
	}

	// loading extra parameters
	pm::DistanceRange range;
	float sigma = 3.0f, penFactor = 3.0f, minAngle = M_PI / 16;
	pm::Lattice::Measure measure = pm::Lattice::FMeasure;
	float beta = 0.1f;
	bool approx = true, stable = false;
	unsigned int N = 100;
	if(options != NULL){
		if(!mxIsStruct(options)){
			mexErrMsgIdAndTxt("MATLAB:args:options", "Options should be a struct!");
		}
		const mxArray *tmp = NULL;
		if(tmp = mxGetField(options, 0, "range")) {
			std::vector<float> R;
			mxLoadVector(R, tmp, "Invalid range vector!");
			if (R.size() >= 1) range.min = R[0];
			if (R.size() >= 2) range.max = R[1];
		}
		if(tmp = mxGetField(options, 0, "sigma"))
			sigma = mxCheckedScalar(tmp, "Invalid sigma parameter!");
		if(tmp = mxGetField(options, 0, "penalty"))
			penFactor = mxCheckedScalar(tmp, "Invalid penalty parameter!");
		if(tmp = mxGetField(options, 0, "min_angle"))
			minAngle = mxCheckedScalar(tmp, "Invalid min_angle parameter!");
		if(tmp = mxGetField(options, 0, "measure")) {
			if(mxStringEquals(tmp, "precision"))
				measure = pm::Lattice::Precision;
			else if(mxStringEquals(tmp, "recall"))
				measure = pm::Lattice::Recall;
			else if(mxStringEquals(tmp, "fmeasure"))
				measure = pm::Lattice::FMeasure;
			else if(mxStringEquals(tmp, "gmeasure"))
				measure = pm::Lattice::GMeasure;
			else if(mxStringEquals(tmp, "value"))
				measure = pm::Lattice::Value;
			else
				mexErrMsgIdAndTxt("MATLAB:lattice:measure", "Invalid measure parameter!");
		}
		if(tmp = mxGetField(options, 0, "beta"))
			beta = mxCheckedScalar(tmp, "Invalid beta parameter!");
		approx = mxBoolField(options, 0, "approximate", true);
		stable = mxBoolField(options, 0, "stable", false);
		if(tmp = mxGetField(options, 0, "N"))
			N = mxCheckedScalar(tmp, "Invalid N parameter!");
	}

	// construct lattice map
	pm::Lattice lattice(points, range, sigma, penFactor);

	switch(mode){
		case Result:
		{	
			pm::GeneratorPair g;
			if(!approx){
				g = lattice.getBestGeneratorPair(minAngle, measure, beta);
			} else {
				g = lattice.approxBestGeneratorPair(N, minAngle, measure, beta);
			}
			
			if(stable){
				pm::Generator g1 = lattice.getStableGenerator(g.p1, minAngle, measure, beta);
				g.p1 = g1.point;
				g.e1 = g1.energy;
				pm::Generator g2 = lattice.getSecondGenerator(g1.point, minAngle, measure, beta);
				g.p2 = g2.point;
				g.e2 = g2.energy;
				
			}

			// store the results
			MatXD G(2, 2, X.type());
			G.update(0, 0, g.p1.x);
			G.update(0, 1, g.p1.y);
			G.update(1, 0, g.p2.x);
			G.update(1, 1, g.p2.y);
			out[0] = G;

			if (nout >= 2) {
				MatXD P(2, 1, IM_MAKETYPE(DataDepth<float>::value, 1));
				P.update(0, 0, g.e1);
				P.update(1, 0, g.e2);
				out[1] = P;
			}
			break;
		}
		case Map:
		{
			// output the mean shift map
			const pm::Bounds &bounds = lattice.boundaries();
			MatXD M(bounds.height(), bounds.width(), IM_MAKETYPE(DataDepth<float>::value, 1));
			for(int y = bounds.min.y; y <= bounds.max.y; ++y) {
				for(int x = bounds.min.x; x <= bounds.max.x; ++x) {
					const float &m = lattice.at(pm::Point2i(x, y));
					M.update(y - bounds.min.y, x - bounds.min.x, m);
				}
			}
			out[0] = M;

			if (nout >= 2) {
				MatXD s(1, 2, IM_MAKETYPE(DataDepth<float>::value, 1));
				s.update(0, 0, bounds.min.x);
				s.update(0, 1, bounds.min.y);
				out[1] = s;
			}
			break;
		}	
		case Energy:
		{
			// output the energy maps
			// from the approximate version
			const pm::Bounds &bounds = lattice.boundaries();
			MatXD E(bounds.height(), bounds.max.x + 1, IM_MAKETYPE(DataDepth<float>::value, 1));
			std::vector<pm::Generator> gens = lattice.getNFirstGenerators(N, measure, beta);
			for(int y = bounds.min.y; y <= bounds.max.y; ++y) {
				for(int x = 0; x <= bounds.max.x; ++x) {
					E.update(y - bounds.min.y, x, 0, 0); // set to zero
				}
			}
			for(int i = 0; i < gens.size(); ++i){
				pm::Generator &g1 = gens[i];
				E.update(g1.point.y - bounds.min.y, g1.point.x, 0, g1.energy);
			}
			out[0] = E;
			
			if(nout > 1){
				MatXD E2(bounds.max.y + 1, bounds.width(), IM_MAKETYPE(DataDepth<float>::value, 1));
				for(int y = 0; y <= bounds.max.y; ++y) {
					for(int x = bounds.min.x; x <= bounds.max.x; ++x) {
						E2.update(y, x - bounds.min.x, 0, 0); // set to zero
					}
				}
				for(int i = 0; i < gens.size(); ++i){
					pm::Generator &g1 = gens[i];
					pm::Generator g2 = lattice.getSecondGenerator(g1.point, minAngle, measure, beta);
					E2.update(g2.point.y, g2.point.x - bounds.min.x, 0, g2.energy);
				}
				out[1] = E2;
			}
			break;
		}
	}
}
