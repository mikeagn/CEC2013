/******************************************************************************
 * Version: 1.0.1
 * Last modified on: 27	October, 2016 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "cec2013.h"

using namespace std;

int main()
{
	/* Demonstration of all functions by 1-D or 2-D plots */
	CEC2013 f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(10);

	/* Create position vectors */
	double x1[1], x2[2];

	/**********************************************************************
	 * Function 1 
	 *********************************************************************/
	cout << "Function 1" << endl;
	fstream file;
	file.open("plots/f1.dat", fstream::out);
	for (x1[0]=0; x1[0]<=30; x1[0]=x1[0]+0.1) {
		file << x1[0] <<"\t"<< f1.evaluate(x1) << endl;
	}
	file.close();

	/**********************************************************************
	 * Function 2 
	 *********************************************************************/
	cout << "Function 2" << endl;
	file.open("plots/f2.dat", fstream::out);
	for (x1[0]=0; x1[0]<=1; x1[0]=x1[0]+0.001) {
		file << x1[0] <<"\t"<< f2.evaluate(x1) << endl;
	}
	file.close();

	/**********************************************************************
	 * Function 3 
	 *********************************************************************/
	cout << "Function 3" << endl;
	file.open("plots/f3.dat", fstream::out);
	for (x1[0]=0; x1[0]<=1; x1[0]=x1[0]+0.001) {
		file << x1[0] <<"\t"<< f3.evaluate(x1) << endl;
	}
	file.close();

	/**********************************************************************
	 * Function 4 
	 *********************************************************************/
	cout << "Function 4" << endl;
	file.open("plots/f4.dat", fstream::out);
	for (x2[0]=-6; x2[0]<=6; x2[0]=x2[0]+0.05) {
		for (x2[1]=-6; x2[1]<=6; x2[1]=x2[1]+0.05) {
			file << f4.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	/**********************************************************************
	 * Function 5 
	 *********************************************************************/
	cout << "Function 5" << endl;
	file.open("plots/f5.dat", fstream::out);
	for (x2[0]=-1.9; x2[0]<=1.9; x2[0]=x2[0]+0.019) {
		for (x2[1]=-1.1; x2[1]<=1.1; x2[1]=x2[1]+0.011) {
			file << f5.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	/**********************************************************************
	 * Function 6 
	 *********************************************************************/
	cout << "Function 6" << endl;
	file.open("plots/f6.dat", fstream::out);
	for (x2[0]=-10; x2[0]<=10; x2[0]=x2[0]+0.1) {
		for (x2[1]=-10; x2[1]<=10; x2[1]=x2[1]+0.1) {
			file << f6.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	/**********************************************************************
	 * Function 7 
	 *********************************************************************/
	cout << "Function 7" << endl;
	file.open("plots/f7.dat", fstream::out);
	for (x2[0]=0.25; x2[0]<=10.001; x2[0]=x2[0]+0.05) {
		for (x2[1]=0.25; x2[1]<=10.001; x2[1]=x2[1]+0.05) {
			file << f7.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	/**********************************************************************
	 * Function 8 
	 *********************************************************************/
	cout << "Function 8" << endl;
	file.open("plots/f8.dat", fstream::out);
	for (x2[0]=0; x2[0]<=1.001; x2[0]=x2[0]+0.01) {
		for (x2[1]=0; x2[1]<=1.001; x2[1]=x2[1]+0.01) {
			file << f8.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	/**********************************************************************
	 * Composition functions F9-F10, 2-D data files for CF1-CF4
	 *********************************************************************/
	cout << "Composition Function 1" << endl;
	int dim = 2;
	CF1 CFunc1(dim);
	file.open("plots/CF1.dat", fstream::out);
	for (x2[0]=-5; x2[0]<=5.001; x2[0]=x2[0]+0.01) {
		for (x2[1]=-5; x2[1]<=5.001; x2[1]=x2[1]+0.01) {
			file << CFunc1.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	cout << "Composition Function 2" << endl;
	CF2 CFunc2(dim);
	file.open("plots/CF2.dat", fstream::out);
	for (x2[0]=-5; x2[0]<=5.001; x2[0]=x2[0]+0.01) {
		for (x2[1]=-5; x2[1]<=5.001; x2[1]=x2[1]+0.01) {
			file << CFunc2.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	cout << "Composition Function 3" << endl;
	CF3 CFunc3(dim);
	file.open("plots/CF3.dat", fstream::out);
	for (x2[0]=-5; x2[0]<=5.001; x2[0]=x2[0]+0.01) {
		for (x2[1]=-5; x2[1]<=5.001; x2[1]=x2[1]+0.01) {
			file << CFunc3.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	cout << "Composition Function 4" << endl;
	CF4 CFunc4(dim);
	file.open("plots/CF4.dat", fstream::out);
	for (x2[0]=-5; x2[0]<=5.001; x2[0]=x2[0]+0.01) {
		for (x2[1]=-5; x2[1]<=5.001; x2[1]=x2[1]+0.01) {
			file << CFunc4.evaluate(x2) <<"\t";
		}
		file << endl;
	}
	file.close();

	return 0;
}
