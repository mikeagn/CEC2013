/******************************************************************************
 * Version: 1.0.1
 * Last modified on: 27	October, 2016 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
#include "cec2013.h"
#include "sortidx.h"

using namespace std;

CEC2013::CEC2013(const int &nfunc) :
	nfunc_(0), cfunc_(NULL), goptima_()
{
	assert(nfunc>0 && nfunc<=20);
	nfunc_ = nfunc;
	init_vars_();

	if (nfunc_ == 11) {
		cfunc_ = new CF1(dimensions_[nfunc_-1]);
	} else if (nfunc_ == 12) {
		cfunc_ = new CF2(dimensions_[nfunc_-1]);
	} else if (nfunc_ == 13 || nfunc_ == 14 || nfunc_ == 16 || nfunc_ == 18 ) {
		cfunc_ = new CF3(dimensions_[nfunc_-1]);
	} else if (nfunc_ == 15 || nfunc_ == 17 || nfunc_ == 19 || nfunc_ == 20 ) {
		cfunc_ = new CF4(dimensions_[nfunc_-1]);
	}

	load_goptima();
}

CEC2013::~CEC2013()
{
	if (cfunc_) delete cfunc_;
}
void CEC2013::init_vars_()
{
	/* specific dimensions for the competition */
	dimensions_[0] = dimensions_[1] = dimensions_[2] = 1;
	dimensions_[3] = dimensions_[4] = dimensions_[5] = dimensions_[6] = 2;
	dimensions_[7] = dimensions_[8] = 3;
	dimensions_[9] = dimensions_[10] = dimensions_[11] = dimensions_[12] = 2;
	dimensions_[13] = dimensions_[14] = 3;
	dimensions_[15] = dimensions_[16] = 5;
	dimensions_[17] = dimensions_[18] = 10;
	dimensions_[19] = 20;
	fopt_[0]  = 200.0;		//1
	fopt_[1]  = 1.0;		//1
	fopt_[2]  = 1.0;		//1
	fopt_[3]  = 200.0;		//2
	fopt_[4]  = 1.031628453489877;	//2
	fopt_[5]  = 186.7309088310239;	//2
	fopt_[6]  = 1.0;		//2
	fopt_[7]  = 2709.093505572820;	//3
	fopt_[8]  = 1.0;		//3
	fopt_[9]  = -2.0;		//2
	for (int i=10; i<20; ++i) { fopt_[i] = 0;}
	rho_[0] = rho_[1] = rho_[2] = rho_[3] = 0.01;
	rho_[4] = rho_[5] = 0.5;
	rho_[6] = 0.2;
	rho_[7] = 0.5;
	rho_[8] = 0.2;
	rho_[9] = rho_[10] = rho_[11] = rho_[12] = rho_[13] = 0.01;
	rho_[14] = rho_[15] = rho_[16] = rho_[17] = rho_[18] = rho_[19] = 0.01;
	nopt_[0] = 2; nopt_[1] = 5; nopt_[2] = 1; nopt_[3] = 4; nopt_[4] = 2;
	nopt_[5] = 18; nopt_[6] = 36; nopt_[7] = 81; nopt_[8] = 216;
	nopt_[9] = 12; nopt_[10] = 6; nopt_[11] = 8; nopt_[12] = 6;
	nopt_[13] = 6; nopt_[14] = 8; nopt_[15] = 6; nopt_[16] = 8;
	nopt_[17] = 6; nopt_[18] = 8; nopt_[19] = 8;

	maxfes_[0] = maxfes_[1] = maxfes_[2] = maxfes_[3] = maxfes_[4] = 50000;
	maxfes_[5] = maxfes_[6] = maxfes_[9] = maxfes_[10] = maxfes_[11] = maxfes_[12] = 200000;
	maxfes_[7] = maxfes_[8] = maxfes_[13] = maxfes_[14] = maxfes_[15] = 400000;
	maxfes_[16] = maxfes_[17] = maxfes_[18] = maxfes_[19] = 400000;
}

double CEC2013::get_lbound(const int &n) const
{
	assert(n>=0 && n<dimensions_[nfunc_-1]);
	double result(0);
	if (nfunc_ == 1 || nfunc_ == 2 || nfunc_ == 3) {
		result = 0;
	} else if (nfunc_ == 4) {
		result = -6;
	} else if (nfunc_ == 5) {
		double tmp[2];
		tmp[0] = -1.9;
		tmp[1] = -1.1;
		result = tmp[n];
	} else if (nfunc_ == 6 || nfunc_ == 8) {
		result = -10;
	} else if (nfunc_ == 7 || nfunc_ == 9) {
		result = 0.25;
	} else if (nfunc_ == 10) {
		result = 0;
	} else if (nfunc_ > 10) {
		result = cfunc_->get_lbound(n);
	}
	return result;
}

double CEC2013::get_ubound(const int &n) const
{
	assert(n>=0 && n<dimensions_[nfunc_-1]);
	double result(0);
	if (nfunc_ == 1) {
		result = 30;
	} else if (nfunc_ == 2 || nfunc_ == 3) {
		result = 1;
	} else if (nfunc_ == 4) {
		result = 6;
	} else if (nfunc_ == 5) {
		double tmp[2];
		tmp[0] = 1.9;
		tmp[1] = 1.1;
		result = tmp[n];
	} else if (nfunc_ == 6 || nfunc_ == 8) {
		result = 10;
	} else if (nfunc_ == 7 || nfunc_ == 9) {
		result = 10;
	} else if (nfunc_ == 10) {
		result = 1;
	} else if (nfunc_ > 10) {
		result = cfunc_->get_ubound(n);
	}
	return result;
}

tFitness CEC2013::evaluate(const double *x)
{
	/* function number: */
	const int dim = dimensions_[nfunc_-1];
	double result(0);
	if (nfunc_ == 1) {
		result =  five_uneven_peak_trap(x, dim);
	} else if (nfunc_ == 2) {
		result =  equal_maxima(x, dim);
	} else if (nfunc_ == 3) {
		result =  uneven_decreasing_maxima(x, dim);
	} else if (nfunc_ == 4) {
		result =  himmelblau(x, dim);
	} else if (nfunc_ == 5) {
		result =  six_hump_camel_back(x, dim);
	} else if (nfunc_ == 6 || nfunc_ == 8) {
		result =  shubert(x, dim);
	} else if (nfunc_ == 7 || nfunc_ == 9) {
		result =  vincent(x, dim);
	} else if (nfunc_ == 10) {
		result =  modified_rastrigin_all(x, dim);
	} else if (nfunc_ > 10) {
		result = cfunc_->evaluate(x);
	}
	return result;
}

tFitness CEC2013::evaluate(std::vector<double> x)
{
	const int dim = x.size();
	assert(dim == dimensions_[nfunc_-1]);
	double *tmpx = new double[dim];
	for (int i=0; i<dim; ++i) { tmpx[i] = x[i]; }

	double result = evaluate(tmpx);
	delete [] tmpx;

	return result;
}

int how_many_goptima(std::vector< std::vector<double> > pop, 
		std::vector< std::vector<double> > &seeds, CEC2013 *pFunc,
		const double &accuracy, const double &radius)
{
	std::vector<std::vector<double> > sorted_pop;
	/* Evaluate pop */
	vector<double> fits;//(pop.size());
	//int i=0;
	for (std::vector< std::vector<double> >::iterator it = pop.begin(); it != pop.end(); ++it) {
		fits.push_back( pFunc->evaluate(*it) );
		//cout << "fits: " << fits[i] << endl;
		//++i;
	}

	/* sorted_pop: sort population by fitness */
	vector<size_t> idx;
	sortIdx(idx, fits, DESCEND);
	for (size_t i=0; i< pop.size(); ++i) {
		sorted_pop.push_back( pop[ idx[i] ] );
	}

	/* find seeds in the temp population */
	get_seeds(sorted_pop, seeds, radius);

	/* Based on the accuracy: check which seeds are global optimizers */
	double count(0);
	for (std::vector< std::vector<double> >::iterator sit = seeds.begin(); sit != seeds.end(); ++sit) {
		/* evaluate seed */
		double seed_fitness = pFunc->evaluate(*sit);

		/* |F_seed - F_goptimum| <= accuracy */
		if ( fabs(seed_fitness - pFunc->get_fitness_goptima()) <= accuracy) {
			++count;
		}
		/* save time */
		if (count == pFunc->get_no_goptima()) break;
	}
	return count;
}

void print_pop(std::vector< std::vector<double> > pop)
{
	for (std::vector< std::vector<double> >::iterator it = pop.begin(); 
			it != pop.end(); ++it) {
		for (std::vector<double>::iterator jt = it->begin();
				jt != it->end(); ++jt) {
			cout << *jt << "\t";
		}
		cout << endl;
	}
}

/**
 * Input:  population vector of vectors,
 * Output: seeds
 */
void get_seeds(std::vector< std::vector<double> > &pop, std::vector< std::vector<double> > &seeds,
		const double &radius)
{
	/* remove any forgoten items in the seed list */
	seeds.clear();
	/* Determine the species seeds: 
	 * iterate through sorted population */
	for (std::vector< std::vector<double> >::iterator it = pop.begin(); it != pop.end(); ++it) {
		bool found(false);
		/* Iterate seeds */
		//cout << "-----" << endl;
		for (std::vector< std::vector<double> >::iterator sit = seeds.begin(); sit != seeds.end(); ++sit) {
			/* Calculate distance from seeds */
			double dist = get_eudist(*it, *sit);
			//cout << "DIST: " << dist << endl;
			/* If the Euclidean distance is less than the radius */
			if (dist <= radius) {
				found = true;
				break;
			}
		} //seeds

		/* If it is not similar to any other seed, then it is a new seed */
		if (!found) { seeds.push_back(*it); }
	} //pop
}

double get_eudist(std::vector<double> &v1, std::vector<double> &v2)
{
	assert(v1.size() == v2.size());
	//if (v1.size() != v2.size()) { exit(-1);}
	double res(0);
	for (size_t i=0; i < v1.size(); i++) {
		res += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	return sqrt(res);
}

void CEC2013::load_goptima() 
{
	if (nfunc_ > 10) {
		goptima_ = cfunc_->get_copy_of_goptima();
		return;
	} 
	std::string filename;
	if (nfunc_ == 1) {
		filename = "data/F1_opt.dat";
	} else if (nfunc_ == 2) {
		filename = "data/F2_opt.dat";
	} else if (nfunc_ == 3) {
		filename = "data/F3_opt.dat";
	} else if (nfunc_ == 4) {
		filename = "data/F4_opt.dat";
	} else if (nfunc_ == 5) {
		filename = "data/F5_opt.dat";
	} else if (nfunc_ == 6) {
		filename = "data/F6_2D_opt.dat";
	} else if (nfunc_ == 8) {
		filename = "data/F6_3D_opt.dat";
	} else if (nfunc_ == 7) {
		filename = "data/F7_2D_opt.dat";
	} else if (nfunc_ == 9) {
		filename = "data/F7_3D_opt.dat";
	} else if (nfunc_ == 10) {
		filename = "data/F8_2D_opt.dat";
	}
	load_goptima(filename);
}

void CEC2013::load_goptima(const std::string &filename)
{
	std::ifstream file;
	file.open(filename.c_str(), std::fstream::in);
	if (!file.is_open()) {
		cerr<< "Error: Can not open file: " << filename << endl;
		//exit(-1);
	}
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		vector<double> indi;
		double tmp;
		for (int i=0; i<dimensions_[nfunc_-1]; ++i) {
			iss >> tmp;
			indi.push_back(tmp);
			//cout << tmp << "\t";
		}
		goptima_.push_back(indi);
		//cout << endl;
	}
	// Close file
	file.close();
}

std::vector< std::vector<double> > CEC2013::get_copy_of_goptima() const
{
	return std::vector< std::vector<double> >(goptima_);
}
