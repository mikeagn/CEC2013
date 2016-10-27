/******************************************************************************
 * Version: 1.0.1
 * Last modified on: 27	October, 2016 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
#include "cfunction.h"

using namespace std;

/******************************************************************************
 * Some declarations
 *****************************************************************************/
/* 
 * return  1      % Minimization 
 * return -1      % Maximization
 * */
static int minmax() { return -1; }
static double E = exp(1);

/******************************************************************************
 * Composition Functions
 *****************************************************************************/
/* Constructors */
CFunction::CFunction() : 
	dimension_(-1), nofunc_(-1), C_(-1), lambda_(NULL), sigma_(NULL),
	bias_(NULL), O_(NULL), M_(NULL), weight_(NULL), lbound_(NULL),
	ubound_(NULL), fi_(NULL), z_(NULL), f_bias_(0), fmaxi_(NULL), 
	tmpx_(NULL), function_(NULL)
{
}

CFunction::CFunction(const int &dim, const int &nofunc) : 
	dimension_(dim), nofunc_(nofunc), C_(2000.0), lambda_(NULL), 
	sigma_(NULL), bias_(NULL), O_(NULL), M_(NULL), weight_(NULL),
	lbound_(NULL), ubound_(NULL), fi_(NULL), z_(NULL), f_bias_(0),
	fmaxi_(NULL), tmpx_(NULL), function_(NULL)
{
	allocate_memory();
}

/* Destructor */
CFunction::~CFunction()
{
	deallocate_memory();
}

void CFunction::allocate_memory()
{
	if (dimension_ < 1 || nofunc_ < 1) return;
	lambda_ = new double[nofunc_];
	sigma_	= new double[nofunc_];
	bias_	= new double[nofunc_];
	O_ 	= new double*[nofunc_];
	M_ 	= new double**[nofunc_];
	for (int i=0; i<nofunc_; ++i) {
		O_[i] = new double[dimension_];
		M_[i] = new double*[dimension_];
		for (int j=0; j<dimension_; ++j) {
			M_[i][j] = new double[dimension_];
		}
	}
	weight_	= new double[nofunc_];
	lbound_	= new double[dimension_];
	ubound_	= new double[dimension_];
	fi_	= new double[nofunc_];
	z_	= new double[dimension_];
	fmaxi_	= new double[nofunc_];
	tmpx_	= new double[dimension_];
	function_ = new compFunction[nofunc_];
}

void CFunction::deallocate_memory()
{
	if (dimension_ < 1 || nofunc_ < 1) return;
	/* Clean up the mess */
	delete [] function_;
	delete [] tmpx_;
	delete [] fmaxi_;
	delete [] z_;
	delete [] fi_;
	delete [] ubound_;
	delete [] lbound_;
	delete [] weight_;
	for (int i=0; i<nofunc_; ++i) { 
		for (int j=0; j<dimension_; ++j) {
			delete [] M_[i][j];
		}
		delete [] M_[i]; 
	}
	delete [] M_;
	for (int i=0; i<nofunc_; ++i) { delete [] O_[i]; }
	delete [] O_;
	delete [] bias_;
	delete [] sigma_;
	delete [] lambda_;
}

void CFunction::calculate_weights(const double *x)
{
	double sum(0), maxi(-INF), maxindex(0);
	for (int i=0; i<nofunc_; ++i) {
		sum = 0.0;
		for (int j=0; j<dimension_; ++j) {
			sum += ( x[j] - O_[i][j] ) * ( x[j] - O_[i][j] );
		}
		weight_[i] = exp( -sum/(2.0 * dimension_ * sigma_[i] * sigma_[i]) );
		if (i==0) { maxi = weight_[i]; }
		if (weight_[i] > maxi) {
			maxi = weight_[i];
			maxindex = i;
		}
		//maxi = max(maxi, weight_[i]);
	}
	sum = 0.0;
	for (int i=0; i<nofunc_; ++i) {
		//if (weight_[i] != maxi) {
		if (i != maxindex) {
			weight_[i] *= (1.0 - pow(maxi, 10.0));
		}
		sum += weight_[i];
	}
	for (int i=0; i<nofunc_; ++i) {
		if (sum == 0.0) {
			weight_[i] = 1.0/(double)nofunc_;
		} else {
			weight_[i] /= sum;
		}
	}
//	for (int i=0; i<nofunc_; ++i) {
//		cout << weight_[i] << "\t";
//	}
//	cout << endl;
}

void CFunction::load_optima(const std::string &filename)
{
	fstream file;
	file.open(filename.c_str(), std::fstream::in);
	if (!file.is_open()) {
		cerr<< "Error: Can not open file: " << filename << endl;
		exit(0);
	}
	double tmp;
	for (int i=0; i< nofunc_; ++i) {
		for (int j=0; j< dimension_; ++j) {
			file >> tmp; 
			O_[i][j] = tmp;
//			cout << O_[i][j] << "\t";
		}
//		cout << endl;
	}
	file.close();
}

void CFunction::load_rotmat(const std::string &filename)
{
	fstream file;
	file.open(filename.c_str(), std::fstream::in);
	if (!file.is_open()) {
		cerr<< "Error: Can not open file: " << filename << endl;
		exit(0);
	}
	double tmp(-1);
	for (int i=0; i<nofunc_; ++i) {
//		cout << "Matrix: " << i << endl;
		for (int j=0; j<dimension_; ++j) {
			for (int k=0; k<dimension_; ++k) {
				file >> tmp; 
				M_[i][j][k] = tmp;
//				cout << M_[i][j][k] << "\t";
			}
//			cout <<endl;
		}
//		cout << endl;
	}
	file.close();
}

void CFunction::init_rotmat_identity()
{
	for (int i=0; i<nofunc_; ++i) {
//		cout << "Matrix: " << i << endl;
		for (int j=0; j<dimension_; ++j) {
			for (int k=0; k<dimension_; ++k) {
				M_[i][j][k] = (j==k ? 1 : 0 );
//				cout << M_[i][j][k] << "\t";
			}
//			cout <<endl;
		}
//		cout << endl;
	}
}

void CFunction::init_optima_rand()
{
	for (int i=0; i< nofunc_; ++i) {
		for (int j=0; j< dimension_; ++j) {
			O_[i][j] = lbound_[j] + (ubound_[j] - lbound_[j]) * rand_uniform();
//			cout << O_[i][j] << "\t";
		}
//		cout << endl;
	}
}

void CFunction::transform_to_z(const double *x, const int &index)
{
	/* Calculate z_i = (x - o_i)/\lambda_i */
	for (int i=0; i<dimension_; ++i) {
		tmpx_[i] = (x[i] - O_[index][i])/lambda_[index];
	}
	/* Multiply z_i * M_i */
	for (int i=0; i<dimension_; ++i) {
		z_[i] = 0;
		for (int j=0; j<dimension_; ++j) {
			/* in MATLAB: M.M1*tmpx' */
			//z_[i] += M_[index][i][j] * tmpx_[j];

			/* in MATLAB: tmpx*M.M1 */
			z_[i] += M_[index][j][i] * tmpx_[j];
		}
//		cout << "i: " << i << " "<< tmpx_[i] << " " << x[i] << " " 
//			<< O_[index][i] <<  " " << lambda_[index] << " "
//			<< z_[i] << endl;
	}
}

void CFunction::transform_to_z_noshift(const double *x, const int &index)
{
	/* Calculate z_i = (x - o_i)/\lambda_i */
	for (int i=0; i<dimension_; ++i) {
		//tmpx_[i] = (x[i] - O_[index][i])/lambda_[index];
		tmpx_[i] = (x[i])/lambda_[index];
	}
	/* Multiply z_i * M_i */
	for (int i=0; i<dimension_; ++i) {
		z_[i] = 0;
		for (int j=0; j<dimension_; ++j) {
			/* in MATLAB: M.M1*tmpx' */
			//z_[i] += M_[index][i][j] * tmpx_[j];

			/* in MATLAB: tmpx*M.M1 */
			z_[i] += M_[index][j][i] * tmpx_[j];
		}
	}
}

void CFunction::calculate_fmaxi()
{
	/* functions */
	for (int i=0; i<nofunc_; ++i) assert(function_[i] != NULL);
	double *x5 = new double[dimension_];
	for (int i=0; i<dimension_; ++i) { x5[i] = 5 ; }

	for (int i=0; i<nofunc_; ++i) {
		transform_to_z_noshift(x5, i);
		fmaxi_[i] = (*function_[i])(z_, dimension_);
//		cout << "FMAXI: " << i << " : " << fmaxi_[i] << endl;
	}
	delete [] x5;
}

tFitness CFunction::evaluate_inner_(const double *x)
{
	double result(0);
	calculate_weights(x);
	for (int i=0; i<nofunc_; ++i) {
		transform_to_z(x, i);
		fi_[i] = (*function_[i])(z_, dimension_);
//		cout << "Func: " << i << " : " << fi_[i] << endl;
	}
	for (int i=0; i<nofunc_; ++i) {
		result += weight_[i]*( C_ * fi_[i] / fmaxi_[i] + bias_[i] );
	}

	return result * minmax() + f_bias_;
}

std::vector< std::vector<double> > CFunction::get_copy_of_goptima() const
{
	assert(O_ != NULL && "O_ == NULL");
	std::vector< std::vector<double> > OO;

	for (int i=0; i< nofunc_; ++i) {
		std::vector<double> kk;
		for (int j=0; j< dimension_; ++j) {
			kk.push_back(O_[i][j]);
		}
		OO.push_back(kk);
	}
	return OO;
}


CF1::CF1(const int dim) : CFunction(dim, 6)
{
	for (int i=0; i<nofunc_; ++i) {
		sigma_[i] = 1;
		bias_[i]  = 0;
		weight_[i]= 0;
	}
	lambda_[0] = 1.0; 
	lambda_[1] = 1.0; 
	lambda_[2] = 8.0; 
	lambda_[3] = 8.0; 
	lambda_[4] = 1.0/5.0; 
	lambda_[5] = 1.0/5.0;
	/* Lower/Upper Bounds */
	for (int i=0; i<dimension_; ++i) {
		lbound_[i] = -5.0;
		ubound_[i] = 5.0;
	}
	/* load optima */
	if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
			|| dimension_ == 10 || dimension_ == 20 ) {
		std::string fname;
		fname = "data/CF1_M_D" + number_to_string(dim) + "_opt.dat";
		load_optima(fname);
//		fname = "data/CF3_M_D" + number_to_string(dim) + ".dat";
//		load_rotmat(fname);
	} else { 
		init_optima_rand();
	}
	/* M_ Identity matrices */
	init_rotmat_identity();
	/* Initialize functions of the composition */
	function_[0] = function_[1] = &FGriewank;
	function_[2] = function_[3] = &FWeierstrass;
	function_[4] = function_[5] = &FSphere;
	//TODO: calculate this correctly in the initialization phase
	calculate_fmaxi();
}

tFitness CF1::evaluate(const double *x)
{
	return evaluate_inner_(x);
}

CF2::CF2(const int dim) : CFunction(dim, 8)
{
	for (int i=0; i<nofunc_; ++i) {
		sigma_[i] = 1.0;
		bias_[i]  = 0.0;
		weight_[i]= 0.0;
	}
	lambda_[0] = 1.0; 
	lambda_[1] = 1.0; 
	lambda_[2] = 10.0; 
	lambda_[3] = 10.0; 
	lambda_[4] = 1.0/10.0; 
	lambda_[5] = 1.0/10.0;
	lambda_[6] = 1.0/7.0;
	lambda_[7] = 1.0/7.0;
	/* Lower/Upper Bounds */
	for (int i=0; i<dimension_; ++i) {
		lbound_[i] = -5.0;
		ubound_[i] = 5.0;
	}
	/* load optima */
	if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
			|| dimension_ == 10 || dimension_ == 20 ) {
		std::string fname;
		fname = "data/CF2_M_D" + number_to_string(dim) + "_opt.dat";
		load_optima(fname);
	} else { 
		init_optima_rand();
	}
	/* M_ Identity matrices */
	init_rotmat_identity();

	/* Initialize functions of the composition */
	function_[0] = function_[1] = &FRastrigin;
	function_[2] = function_[3] = &FWeierstrass;
	function_[4] = function_[5] = &FGriewank;
	function_[6] = function_[7] = &FSphere;
	calculate_fmaxi();
}

tFitness CF2::evaluate(const double *x)
{
	return evaluate_inner_(x);
}

CF3::CF3(const int dim) : CFunction(dim, 6)
{
	for (int i=0; i<nofunc_; ++i) {
		bias_[i]  = 0.0;
		weight_[i]= 0.0;
	}
	sigma_[0] = 1.0;
	sigma_[1] = 1.0;
	sigma_[2] = 2.0;
	sigma_[3] = 2.0;
	sigma_[4] = 2.0;
	sigma_[5] = 2.0;
	lambda_[0] = 1.0/4.0; 
	lambda_[1] = 1.0/10.0; 
	lambda_[2] = 2.0; 
	lambda_[3] = 1.0; 
	lambda_[4] = 2.0; 
	lambda_[5] = 5.0;
	/* Lower/Upper Bounds */
	for (int i=0; i<dimension_; ++i) {
		lbound_[i] = -5.0;
		ubound_[i] = 5.0;
	}
	/* load optima */
	if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
			|| dimension_ == 10 || dimension_ == 20 ) {
		std::string fname;
		fname = "data/CF3_M_D" + number_to_string(dim) + "_opt.dat";
		load_optima(fname);
		fname = "data/CF3_M_D" + number_to_string(dim) + ".dat";
		load_rotmat(fname);
	} else { 
		init_optima_rand();
		//TODO: Generate dimension independent rotation matrices
		/* M_ Identity matrices */
		init_rotmat_identity();
	}
	/* Initialize functions of the composition */
	function_[0] = function_[1] = &FEF8F2;
	function_[2] = function_[3] = &FWeierstrass;
	function_[4] = function_[5] = &FGriewank;
	calculate_fmaxi();
}

tFitness CF3::evaluate(const double *x)
{
	return evaluate_inner_(x);
}

CF4::CF4(const int dim) : CFunction(dim, 8)
{
	for (int i=0; i<nofunc_; ++i) {
		sigma_[i] = 1.0;
		bias_[i]  = 0.0;
		weight_[i]= 0.0;
	}
	sigma_[0] = 1.0;
	sigma_[1] = 1.0;
	sigma_[2] = 1.0;
	sigma_[3] = 1.0;
	sigma_[4] = 1.0;
	sigma_[5] = 2.0;
	sigma_[6] = 2.0;
	sigma_[7] = 2.0;
	lambda_[0] = 4.0; 
	lambda_[1] = 1.0; 
	lambda_[2] = 4.0; 
	lambda_[3] = 1.0; 
	lambda_[4] = 1.0/10.0; 
	lambda_[5] = 1.0/5.0;
	lambda_[6] = 1.0/10.0;
	lambda_[7] = 1.0/40.0;
	/* Lower/Upper Bounds */
	for (int i=0; i<dimension_; ++i) {
		lbound_[i] = -5.0;
		ubound_[i] = 5.0;
	}
	/* load optima */
	if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
			|| dimension_ == 10 || dimension_ == 20) {
		std::string fname;
		fname = "data/CF4_M_D" + number_to_string(dim) + "_opt.dat";
		load_optima(fname);
		fname = "data/CF4_M_D" + number_to_string(dim) + ".dat";
		load_rotmat(fname);
	} else {
		init_optima_rand();
		//TODO: make dimension independent rotation matrices
		/* M_ Identity matrices */
		init_rotmat_identity();
	}
	/* Initialize functions of the composition */
	function_[0] = function_[1] = &FRastrigin;
	function_[2] = function_[3] = &FEF8F2;
	function_[4] = function_[5] = &FWeierstrass;
	function_[6] = function_[7] = &FGriewank;
	calculate_fmaxi();
}

tFitness CF4::evaluate(const double *x)
{
	return evaluate_inner_(x);
}

/******************************************************************************
 * Basic Benchmark functions 
 *****************************************************************************/
/******************************************************************************
 * F1: Five-Uneven-Peak Trap 
 * Variable ranges: x in [0, 30
 * No. of global peaks: 2
 * No. of local peaks:  3. 
 *****************************************************************************/
tFitness five_uneven_peak_trap(const double *x, const int &dim)
{
	tFitness result=-1.0;
	if (x[0]>=0 && x[0]< 2.5) {
		result = 80*(2.5-x[0]);
	} else if (x[0] >= 2.5 && x[0] < 5.0) {
		result = 64*(x[0]-2.5);
	} else if (x[0] >= 5.0 && x[0] < 7.5) {
		result = 64*(7.5-x[0]);
	} else if (x[0] >= 7.5 && x[0] < 12.5) {
		result = 28*(x[0]-7.5);
	} else if (x[0] >= 12.5 && x[0] < 17.5) {
		result = 28*(17.5-x[0]);
	} else if (x[0] >= 17.5 && x[0] < 22.5) {
		result = 32*(x[0]-17.5);
	} else if (x[0] >= 22.5 && x[0] < 27.5) {
		result = 32*(27.5-x[0]);
	} else if (x[0] >= 27.5 && x[0] <= 30) {
		result = 80*(x[0]-27.5);
	}

	return result;
}

/******************************************************************************
 * F2: Equal Maxima
 * Variable ranges: x in [0, 1]
 * No. of global peaks: 5
 * No. of local peaks:  0. 
 *****************************************************************************/
tFitness equal_maxima(const double *x, const int &dim)
{
	tFitness s = sin(5.0 * M_PI * x[0]);
	return pow(s, 6);
}

/******************************************************************************
 * F3: Uneven Decreasing Maxima
 * Variable ranges: x in [0, 1]
 * No. of global peaks: 1
 * No. of local peaks:  4. 
 *****************************************************************************/
tFitness uneven_decreasing_maxima(const double *x, const int &dim)
{
	tFitness tmp1 = -2*log(2)*((x[0]-0.08)/0.854)*((x[0]-0.08)/0.854);
	tFitness tmp2 = sin( 5*M_PI*(pow(x[0],3.0/4.0)-0.05) );
	return exp(tmp1) * pow(tmp2, 6);
}

/******************************************************************************
 * F4: Himmelblau
 * Variable ranges: x, y in [−6, 6
 * No. of global peaks: 4
 * No. of local peaks:  0.
 *****************************************************************************/
tFitness himmelblau(const double *x, const int &dim)
{
	return 200 - (x[0]*x[0] + x[1] - 11)*(x[0]*x[0] + x[1] - 11) - 
		(x[0] + x[1]*x[1] - 7)*(x[0] + x[1]*x[1] - 7);
}

/******************************************************************************
 * F5: Six-Hump Camel Back
 * Variable ranges: x in [−1.9, 1.9]; y in [−1.1, 1.1]
 * No. of global peaks: 2
 * No. of local peaks:  2.
 *****************************************************************************/
tFitness six_hump_camel_back(const double *x, const int &dim)
{
	return -( (4 - 2.1*x[0]*x[0] + pow(x[0],4.0)/3.0)*x[0]*x[0] + 
		x[0]*x[1] + (4*x[1]*x[1] -4)*x[1]*x[1] );
}

/******************************************************************************
 * F6: Shubert
 * Variable ranges: x_i in  [−10, 10]^n, i=1,2,...,n
 * No. of global peaks: n*3^n
 * No. of local peaks: many
 *****************************************************************************/
tFitness shubert(const double *x, const int &dim)
{
	tFitness result(1), sum(0); 
	for (int i=0; i<dim; i++) {
		sum=0;
		for (int j=1; j<6; j++) {
			sum = sum + j * cos((j+1) * x[i] + j);
		}
		result = result * sum;
	}
	return -result;
}

/******************************************************************************
 * F7: Vincent
 * Variable range: x_i in [0.25, 10]^n, i=1,2,...,n
 * No. of global optima: 6^n
 * No. of local optima:  0.
 *****************************************************************************/
tFitness vincent(const double *x, const int &dim)
{
	tFitness result(0);
	for (int i=0; i<dim; i++){
		if (x[i]<=0){
			cerr << "Illegal value: " << x[i] << endl;
			exit(-1);
		}
		result = result + sin(10 * log(x[i]));
	}
	return result/dim;
}

/******************************************************************************
 * F8: Modified Rastrigin - All Global Optima
 * Variable ranges: x_i in [0, 1]^n, i=1,2,...,n
 * No. of global peaks: \prod_{i=1}^n k_i
 * No. of local peaks:  0.
 *****************************************************************************/
/* Modified Rastrigin -- All Global Optima */
static double MPPF92[2] = {3, 4};
static double MPPF98[8] = {1, 2, 1, 2, 1, 3, 1, 4};
static double MPPF916[16] = {1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 4};

tFitness modified_rastrigin_all(const double *x, const int &dim)
{
	tFitness result(0);
	for (int i=0; i<dim; i++){
		if (dim == 2)  { result = result + 10+ 9*cos(2*M_PI*MPPF92[i]*x[i]); }
		if (dim == 8)  { result = result + 10+ 9*cos(2*M_PI*MPPF98[i]*x[i]); }
		if (dim == 16) { result = result + 10+ 9*cos(2*M_PI*MPPF916[i]*x[i]); }
	}
	return -result;
}

/******************************************************************************
 * Basic functions for composition 
 *****************************************************************************/
/* Ackley's function */
tFitness FAckley(const double *x, const int &dim)
{
	tFitness sum1(0.0), sum2(0.0), result;
	for (int i=0; i<dim; ++i) {
		sum1 += x[i]*x[i];
		sum2 += cos(2.0*M_PI*x[i]);
	}
	sum1 = -0.2*sqrt(sum1/dim);
	sum2 /= dim;
	result = 20.0 + E - 20.0*exp(sum1) - exp(sum2);
	return result;
}

/* Rastrigin's function */
tFitness FRastrigin(const double *x, const int &dim)
{
	tFitness result(0.0);
	for (int i=0; i<dim; ++i) {
		result += (x[i]*x[i] - 10.0*cos(2.0*M_PI*x[i]) + 10.0);
	}
	return result;
}

/* Weierstrass's function */
tFitness FWeierstrass(const double *x, const int &dim)
{
	tFitness result(0.0), sum(0.0), sum2(0.0), a(0.5), b(3.0);
	int k_max(20);

	for (int j=0; j<=k_max; ++j) {
		sum2 += pow(a,j)*cos(2.0*M_PI*pow(b,j)*(0.5));
	}
	for (int i=0; i<dim; ++i) {
		sum = 0.0;
		for (int j=0; j<=k_max; ++j) {
			sum += pow(a,j)*cos(2.0*M_PI*pow(b,j)*(x[i]+0.5));
		}
		result += sum;
	}
	return result - sum2*dim;
}

/* Griewank's function */
tFitness FGriewank(const double *x, const int &dim)
{
	tFitness sum(0.0), prod(1.0), result(0.0);

	for (int i=0; i<dim; ++i) {
		sum  += x[i]*x[i]/4000.0;
		prod *= cos( x[i]/sqrt(double(1.0+i)) );
	}
	result = 1.0 + sum - prod;
	return result;
}

/* Sphere function */
tFitness FSphere(const double *x, const int &dim)
{
	tFitness result(0.0);
	for (int i=0; i<dim; ++i) {
		result += x[i]*x[i];
	}
	return result;
}

/* Schwefel's function */
tFitness FSchwefel(const double *x, const int &dim)
{
	tFitness sum1(0.0), sum2(0.0);

	for (int i=0; i<dim; ++i) {
		sum2 = 0.0;
		for (int j=0; j<=i; ++j) {
			sum2 += x[j];
		}
		sum1 += sum2*sum2;
	}
	return sum1;
}

/* Rosenbrock's function */
tFitness FRosenbrock(const double *x, const int &dim)
{
	tFitness result(0.0);

	for (int i=0; i<dim-1; ++i) {
		result += 100.0*pow((x[i]*x[i]-x[i+1]),2.0) + 1.0*pow((x[i]-1.0),2.0);
	}
	return result;
}

/* FEF8F2 function */
tFitness FEF8F2(const double *xx, const int &dim)
{
	tFitness result(0.0);
	double x(0), y(0), f(0), f2(0);

	for (int i=0; i<dim-1; ++i) {
		x = xx[i]   +1;
		y = xx[i+1] +1;

		f2 = 100.0*(x*x - y)*(x*x - y) + (1.0 - x)*(1.0 - x);
		f  = 1.0 + f2*f2/4000.0 - cos(f2);

		result += f;
	}
	/* do not forget the (dim-1,0) case! */
	x = xx[dim-1] +1;
	y = xx[0]     +1;

	f2 = 100.0*(x*x - y)*(x*x - y) + (1.0 - x)*(1.0 - x);
	f  = 1.0 + f2*f2/4000.0 - cos(f2);

	result += f;

	return result;
}

