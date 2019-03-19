/******************************************************************************
 * Version: 1.0.1
 * Last modified on: 27	October, 2016 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
#ifndef __CEC2013_FUNCTIONS_H__
#define __CEC2013_FUNCTIONS_H__
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
//#include <values.h>
#include <limits.h>
#include <float.h>
#include <cassert>
//TODO: random number generator
#include "rand2.h"

/* Just to define M_PI */
#define _USE_MATH_DEFINES
//#define INF MAXDOUBLE
#define INF DBL_MAX

/* Define tFitness */
#ifndef _TFITNESS 
#define _TFITNESS 1
// This should be a real number i.e. float, double, long double
typedef long double tFitness; 	// Fitness type
//typedef long double varType; 	// Variable type
#endif

/* type definition for an array of function pointers */
typedef tFitness (*compFunction)(const double *, const int &);

/* Composition Functions framework */
class CFunction
{
	//non-copyable
	CFunction(const CFunction&);
	CFunction& operator=(const CFunction &);
public:
	CFunction();
	CFunction(const int &dim, const int &nofunc);
	virtual ~CFunction();

	virtual tFitness evaluate(const double *x) = 0;
	double get_lbound(const int &ivar) const { return lbound_[ivar]; } 
	double get_ubound(const int &ivar) const { return ubound_[ivar]; } 

protected:
//public:
	int dimension_;
	int nofunc_;
	double C_;// = 2000.0;
	double *lambda_;
	double *sigma_;
	double *bias_;
	double **O_;
	double ***M_;
	double *weight_;
	double *lbound_;
	double *ubound_;
	double *fi_;
	double *z_;
	double f_bias_;
	double *fmaxi_;
	double *tmpx_;
	compFunction *function_;

	/* Inner help functions */
	void allocate_memory();
	void deallocate_memory();
	void init_rotmat_identity();
	void init_optima_rand();
	void load_optima(const std::string &filename);
	void load_rotmat(const std::string &filename);
	void calculate_weights(const double *x);
	void transform_to_z(const double *x, const int &index);
	void transform_to_z_noshift(const double *x, const int &index);
	void calculate_fmaxi();
	tFitness evaluate_inner_(const double *x);
	std::vector< std::vector<double> > get_copy_of_goptima() const;
	friend class CEC2013;
};
/* Basic Benchmark functions */
tFitness five_uneven_peak_trap(const double *x, const int &dim);
tFitness equal_maxima(const double *x, const int &dim);
tFitness uneven_decreasing_maxima(const double *x, const int &dim);
tFitness himmelblau(const double *x, const int &dim);
tFitness six_hump_camel_back(const double *x, const int &dim);
tFitness shubert(const double *x, const int &dim);
tFitness vincent(const double *x, const int &dim);
tFitness modified_rastrigin_all(const double *x, const int &dim);

/* Basic functions for composition */
tFitness FSphere(const double *x, const int &dim);
tFitness FAckley(const double *x, const int &dim);
tFitness FEF8F2(const double *xx, const int &dim);
tFitness FGriewank(const double *x, const int &dim);
tFitness FSchwefel(const double *x, const int &dim);
tFitness FRastrigin(const double *x, const int &dim);
tFitness FRosenbrock(const double *x, const int &dim);
tFitness FWeierstrass(const double *x, const int &dim);

/* Interfaces for Composition functions */
class CF1 : public CFunction
{
	//non-copyable
	CF1(const CF1 &);
	CF1& operator=(const CF1&);
public:
	CF1(const int dim);
	tFitness evaluate(const double *x);
};

class CF2 : public CFunction
{
	//non-copyable
	CF2(const CF2 &);
	CF2& operator=(const CF2&);
public:
	CF2(const int dim);
	tFitness evaluate(const double *x);
};

class CF3 : public CFunction
{
	//non-copyable
	CF3(const CF3 &);
	CF3& operator=(const CF3&);
public:
	CF3(const int dim);
	tFitness evaluate(const double *x);
};

class CF4 : public CFunction
{
	//non-copyable
	CF4(const CF4 &);
	CF4& operator=(const CF4&);
public:
	CF4(const int dim);
	tFitness evaluate(const double *x);
};

/* Help template to fix file names */
template <typename T>
std::string number_to_string(T no)
{
	std::ostringstream ss; ss << no; return ss.str();
}
#endif
