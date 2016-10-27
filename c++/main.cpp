/******************************************************************************
 * Version: 1.0.1
 * Last modified on: 27	October, 2016 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
#include <iostream>
#include <cstdlib>
#include <vector>

#include "cec2013.h"

using namespace std;

int main()
{
	/* Demonstration of all functions */
	CEC2013 f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8), 
		f9(9), f10(10), f11(11), f12(12), f13(13), f14(14), 
		f15(15), f16(16), f17(17), f18(18), f19(19), f20(20);

	/* Create position vectors */
	double x1[1] = {1}; 
	double x2[2] = {1,1};
	double x3[3] = {1,1,1};
	double x5[5] = {1,1,1,1,1};
	double x10[10];
	double x20[20];
	for (int i=0; i<10; i++) { x10[i] = 1; }
	for (int i=0; i<20; i++) { x20[i] = 1; }

	/* Evaluate :-) */
	cout << "f1(1)      = " << f1.evaluate(x1) << endl;
	cout << "f2(1)      = " << f2.evaluate(x1) << endl;
	cout << "f3(1)      = " << f3.evaluate(x1) << endl;
	cout << "f4(1,1)    = " << f4.evaluate(x2) << endl;
	cout << "f5(1,1)    = " << f5.evaluate(x2) << endl;
	cout << "f6(1,1)    = " << f6.evaluate(x2) << endl;
	cout << "f7(1,1)    = " << f7.evaluate(x2) << endl;
	cout << "f8(1,1,1)  = " << f8.evaluate(x3) << endl;
	cout << "f9(1,1,1)  = " << f9.evaluate(x3) << endl;
	cout << "f10(1,1)   = " << f10.evaluate(x2) << endl;
	cout << "f11(1,1)   = " << f11.evaluate(x2) << endl;
	cout << "f12(1,1)   = " << f12.evaluate(x2) << endl;
	cout << "f13(1,1)   = " << f13.evaluate(x2) << endl;
	cout << "f14(1,1,1) = " << f14.evaluate(x3) << endl;
	cout << "f15(1,1,1) = " << f15.evaluate(x3) << endl;
	cout << "f16(x5)    = " << f16.evaluate(x5) << endl;
	cout << "f17(x5)    = " << f17.evaluate(x5) << endl;
	cout << "f18(x10)   = " << f18.evaluate(x10) << endl;
	cout << "f19(x10)   = " << f19.evaluate(x10) << endl;
	cout << "f20(x20)   = " << f20.evaluate(x20) << endl;

	/* Iterate through functions */
	CEC2013 *pFunc;
	for (int index=1; index<=20; ++index) {
		/* Create one */
		pFunc = new CEC2013(index);

		int dim = pFunc->get_dimension();
		double *X = new double[dim];
		for (int i=0; i<dim; ++i) { X[i] = 1.0; }

		cout << "f" << index << " dim: " << pFunc->get_dimension() 
			<< " f" << index << "(X) = " << pFunc->evaluate(X) << endl;

		/* Clean up */
		delete [] X;
		delete pFunc;
	}

	/**********************************************************************
	 *  Demonstration of using how_many_goptima function 
	 *********************************************************************/
	pFunc = new CEC2013(14);
	/* Create a population: std::vector< std::vector<double> > */
	std::vector< std::vector<double> > pop;
	const int dim(pFunc->get_dimension());
	for (int item=0; item<10; ++item) {
		vector<double> x(dim);
		for (int i=0;i<dim; ++i) {
			x[i] = rand_uniform();
		}
		pop.push_back(x);
	}
	/* Print population */
	cout << "-------------------------------------------------------" <<endl;
	for (std::vector< std::vector<double> >::iterator it = pop.begin(); 
			it != pop.end(); ++it) {
		cout << "Fitness: " << pFunc->evaluate(*it) << "\tGene:\t";
		for (std::vector<double>::iterator jt = it->begin();
				jt != it->end(); ++jt) {
			cout << *jt << "\t";
		}
		cout << endl;
	}
	cout << "-------------------------------------------------------" <<endl;

	/* Calculate how many global optima are in the population */
	double accuracy=0.001;
	std::vector< std::vector<double> > seeds; 
	cout << "In the current population there exist " 
		<< how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) 
		<< " global optimizers." << endl;

//	/* Example usage for getting the global optima *only* for statistical purposes */
//	typedef std::vector< std::vector<double> > Array2D;
//
//	/* Iterate through functions */
//	cout << "Global optima per problem: " << endl;
//	CEC2013 *pFunc_gopt;
//	for (int index=1; index<=20; ++index) {
//		/* Create problem */
//		pFunc_gopt = new CEC2013(index);
//		cout << "Problem: " << index << endl;
//
//		/* Get optima */
//		Array2D goptima = pFunc_gopt->get_copy_of_goptima();
//		cout << "SIZE: " << goptima.size() << endl;
//		for (Array2D::iterator it = goptima.begin(); it != goptima.end(); ++it) {
//			//evaluate optimum
//			cout << "Fitness: " << pFunc_gopt->evaluate(*it) << " f: " << pFunc_gopt->get_fitness_goptima() << "\tGene:\t";
//			for (std::vector<double>::iterator jt = it->begin();
//					jt != it->end(); ++jt) {
//				cout << *jt << "\t";
//			}
//			cout << endl;
//		}
//
//		/* Clean up */
//		delete pFunc_gopt;
//	}

	/* Clean up */
	delete pFunc;
	
	return EXIT_SUCCESS;
}
