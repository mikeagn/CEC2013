/******************************************************************************
 * Version: 1.0
 * Last modified on: 25 April, 2013 
 * Developers: Michael G. Epitropakis
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 * ***************************************************************************/
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "cec2013.h"
#include "params.h"

typedef std::vector< std::vector< std::vector<double> > > Array3D;
typedef std::vector< std::vector<double> > Array2D;

using namespace std;

std::string operator+(const string& a, const string& b);
void load_pop(const std::string filename, std::vector< std::vector<double> > &pop, const int &dim);
void init_memory(Array3D &results, Array2D &PR, Array2D &SR, 
		const int &problems, const int &noruns, const int &noacc);
void gather_stats(Array3D &results, Array2D &PR, Array2D &SR,
		const int &problems, const int &noruns, const int &noacc);
void create_latex_table(Array2D &PR, Array2D &SR, const std::string &outfilename, 
		double *acclevels, const std::string &algname, const int &problems, const int &noacc);

int main(int argc, char **argv)
{
	parameters p(argc, argv);
	const int noruns = p.noruns, problems = p.problems, noacc = 5;
	double acclevels[noacc] = {0.1, 0.01, 0.001, 0.0001, 0.00001};
	std::string outfilename(p.outfilename), respath(p.respath), algname(p.algname);
	Array3D results;
	Array2D PR, SR;
	init_memory(results, PR, SR, problems, noruns, noacc);

	/* Iterate through functions */
	CEC2013 *pFunc;
	for (int index = 1; index <= problems; ++index) {
		/* Create one */
		pFunc = new CEC2013(index);
		const int dim = pFunc->get_dimension();
		cout << "F_" << index << ": ";

		/* for each run */
		for (int irun = 0; irun < noruns; ++irun) {
			cout << ".";
			/* Create and read population from files */
			std::vector< std::vector<double> > pop;
			string tmpf = respath + "problem"  + number_to_string(index) 
								+ "run" + number_to_string(irun+1) + ".dat";

			//problem19.run.50.acc.0.1000000000.dat
		//	string tmpf = respath + "problem"  + number_to_string(index) 
		//						+ ".run." + number_to_string(irun+1) + ".acc.0.0000100000.dat";

			//multi.population.1.s7.acc.0.1000000000.run1000.res
//			string tmpf = respath + "multi.population."  + number_to_string(index) 
//								+ ".s7.acc.0.0000100000.run" + number_to_string(1000+irun) + ".res";
			//multi.population.1.s7.acc.0.1000000000.run1000.res
				//cout << tmpf << endl;
			/* Load population */
			load_pop(tmpf, pop, dim);

			int size = pop.size();
			/* for each accuracy level */
			for (int iacc = 0; iacc < noacc; ++iacc) {
				/* Calculate how many global optima are in the population */
				double accuracy = acclevels[iacc];
				std::vector< std::vector<double> > seeds; 
				if (size == 0) {
					results[index-1][irun][iacc] = 0;
				} else {
					int count = how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho());
					results[index-1][irun][iacc] = count;
				}
			}
		}
		cout << endl;
		/* Clean up */
		delete pFunc;
	}
	gather_stats(results, PR, SR, problems, noruns, noacc);

	/* Create the PR/SR latex table */
	create_latex_table(PR, SR, outfilename, acclevels, algname, problems, noacc);

	return EXIT_SUCCESS;
}

void create_latex_table(Array2D &PR, Array2D &SR, const std::string &outfilename, 
		double *acclevels, const std::string &algname, const int &problems, const int &noacc)
{
	const int rowblocks = 4;
	const int algs_per_rowblock = 5;
	std::fstream file;
	file.open(outfilename.c_str(), std::fstream::out);
	int dims[problems];

	for (int i=0; i<problems; ++i) {
		CEC2013 *pF = new CEC2013(i+1);
		dims[i] = pF->get_dimension();
		delete pF;
	}	

	/* Table's header */
	file << endl;
	file << "\\begin{table}[!ht]" << endl;
	file << "\\centering" << endl;
	file << "\\caption{Peak ratios and success rates of " << algname << " algorithm.}\\label{tab:" << algname << "}" << endl;
	file << "\\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c}" << endl;
	file << "\\hline\\hline" << endl;

	/* Table's body */
	for (int i=0; i<rowblocks; ++i) {
		/* First header */
		file << "Accuracy ";
		for (int col=0; col<algs_per_rowblock; ++col) {
			int prob = (i)*algs_per_rowblock + col ;
			file << "& \\multicolumn{2}{c|}{$F_{"<<prob+1<<"}$ ("<<dims[prob]<<"D)} ";
		}
		file << "\\\\ \\hline" << endl;
		/* Second header */
		file << "level $\\varepsilon$ ";
		for (int col=0; col<algs_per_rowblock; ++col) {
			file << "& PR & SR ";
		}
		file << "\\\\ \\hline" << endl;
		/* for each accuracy level */
		for (int iacc = 0; iacc < noacc; ++iacc) {
			file << std::scientific << std::setprecision(1) << acclevels[iacc] << " ";
			file << std::fixed << std::setprecision(3);
			for (int col=0; col<algs_per_rowblock; ++col) {
				int prob = (i)*algs_per_rowblock + col ;
				file << "& " << PR[prob][iacc] << "\t& " << SR[prob][iacc] << "\t";
			}
			file << "\\\\ " << endl;
		}
		file << " \\hline" << endl;
	}//rowblocks

	/* Table's footer */
	file << "\\hline" << endl;
	file << "\\end{tabular}" << endl;
	file << "\\end{table}" << endl;

	// Close file
	file.close();
}

void init_memory(Array3D &results, Array2D &PR, Array2D &SR, 
		const int &problems, const int &noruns, const int &noacc)
{
	results.resize(problems);
	PR.resize(problems);
	SR.resize(problems);
	for (int i=0; i<problems; ++i) {
		results[i].resize(noruns);
		for (int j=0; j<noruns; ++j) {
			results[i][j].resize(noacc);
		}
		PR[i].resize(noacc);
		SR[i].resize(noacc);
	}
}

void load_pop(const std::string filename, std::vector< std::vector<double> > &pop, const int &dim)
{
	std::ifstream file;
	file.open(filename.c_str(), std::fstream::in);
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		vector<double> indi;
		double tmp;
		for (int i=0; i<dim; ++i) {
			iss >> tmp;
			indi.push_back(tmp);
			//cout << tmp << "\t";
		}
		//cout << endl;
		pop.push_back(indi);
	}
	// Close file
	file.close();
}

void gather_stats(Array3D &results, Array2D &PR, Array2D &SR,
		const int &problems, const int &noruns, const int &noacc)
{
	CEC2013 *pFunc;
	/* for each problem */
	for (int index = 0; index < problems; ++index) {
		pFunc = new CEC2013(index+1);

		/* for each accuracy level */
		for (int iacc = 0; iacc < noacc; ++iacc) {
			double sum1(0); int sum2(0);
			/* for each run */
			for (int irun = 0; irun < noruns; ++irun) {
				sum1 += results[index][irun][iacc];
				if (results[index][irun][iacc] == pFunc->get_no_goptima()) ++sum2;
			}//runs
			PR[index][iacc] = sum1/(noruns*pFunc->get_no_goptima());
			SR[index][iacc] = (1.0*sum2)/(noruns);
		}//acc
		delete pFunc;
	}//problems
}

std::string operator+ (const string& a, const string& b) {
	std::string result(a);
	result += b;
	return result;
}
