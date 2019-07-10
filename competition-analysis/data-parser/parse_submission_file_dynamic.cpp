/******************************************************************************
 * Version: 1.0.2
 * Last modified on: 10 July, 2019 
 * Developers: Michael G. Epitropakis
 *      email: m_(DOT)_epitropakis_(AT)_gmail_(DOT)_com
 * ***************************************************************************/
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <map>
#include <string>

#include "cec2013.h"
#include "params.h"

typedef std::vector< std::vector< std::vector<double> > > Array3D;
typedef std::vector< std::vector<double> > Array2D;
typedef std::pair< std::vector<double>, std::vector<double> > TwoVectors;
typedef std::multimap< std::string, TwoVectors > archive;

using namespace std;

std::string vector_to_string(const std::vector<double> & vec);
std::string operator+(const string& a, const string& b);
void load_file_with_archive_action(const std::string filename, std::vector< std::vector<double> > &data, const int &dim);
int read_data_line(const std::vector<double> &data, TwoVectors &result, const int &dim) ;
void init_memory(Array3D &results, Array2D &PR, Array2D &SR, 
		const int &problems, const int &noruns, const int &noacc);
void calculate_number_of_goptima_per_acc_run_get_result(std::vector<double> &results, const Array2D &pop, 
		const int &noacc, const double *acclevels, CEC2013 *pFunc);

std::string int_to_string_with_leading_zeros(int no, size_t orders_of_magnitute)
{
	std::ostringstream ss; ss << setw(orders_of_magnitute) << setfill('0') << no; return ss.str();
}
static bool header = true;

int main(int argc, char **argv)
{
	parameters p(argc, argv);
	const int noruns = p.noruns, problems = p.problems, noacc = 5;
	double acclevels[noacc] = {0.1, 0.01, 0.001, 0.0001, 0.00001};
	std::string outfilename(p.outfilename), respath(p.respath+"/"), algname(p.algname);
	Array3D results;
	Array2D PR, SR, pop, costs, all_dists, data;
	std::vector<double> fitnesses;
	init_memory(results, PR, SR, problems, noruns, noacc);

	/* for each functions */
	CEC2013 *pFunc;
	for (int index = 1; index <= 20; ++index) {
		/* Create one */
		pFunc = new CEC2013(index);
		const unsigned int dim = pFunc->get_dimension();
		cout << "F_" << index << ": ";
		// load global optima
		Array2D goptima = pFunc->get_copy_of_goptima();

		/* for each run */
		for (int irun = 0; irun < noruns; ++irun) {
			cout << ".";
			/* Create data structures and read data from files */
			pop.clear(); costs.clear(); data.clear();
			//Parse: respath = algorithms/cma/ + problem001run001.dat
			string tmpfilename = respath + "problem"  + int_to_string_with_leading_zeros(index,3) 
								+ "run" + int_to_string_with_leading_zeros(irun+1,3) + ".dat";

			/* Load file */
			load_file_with_archive_action(tmpfilename, data, dim);

			/* sort based on FEs */
			std::sort(data.begin(), data.end(),
					[](const std::vector<double>& a, const std::vector<double>& b) {
							return a[a.size()-3] < b[b.size()-3];
					});

			// create archive
			archive myarchive;
			for (unsigned int ii=0; ii< data.size(); ii++) {
				myarchive.clear();
				for (unsigned int jj=0; jj<=ii; jj++) {
					//cout << data[ii][jj] << "\t";
					TwoVectors result;
					int action = read_data_line(data[jj], result, dim);
					if (action == 0) {
						// empty the archive and 
						myarchive.clear();
						/// add the new solution
						std::string key(vector_to_string(result.first));
						myarchive.insert( std::pair<std::string, TwoVectors>(key, result) );
					} else if (action == 1) {
						// add the new solution
						std::string key(vector_to_string(result.first));
						myarchive.insert( std::pair<std::string, TwoVectors>(key, result) );
					} else if (action == -1) {
						// search for it and remove it
						std::string key(vector_to_string(result.first));
						archive::iterator it = myarchive.find(key);
						if (it != myarchive.end()) { // found it
							myarchive.erase(it);
						} else {
							cerr << "WARNING: I could not remove solution at line " << jj << " because it was not in the archive." << endl;
						}
					} else {
						cerr << "ERROR: " << tmpfilename << " contains incorrect archive action index: " << action << " at line: "<< jj << endl;
						exit(-1);
					}
				}
				Array2D cpop; 
				// add info from archive into pop/costs
				for (archive::iterator it = myarchive.begin(); it != myarchive.end(); ++it) {
					cpop.push_back((it->second).first);
				}

				std::vector<double> res;
				// Calculate how many global optima are in the population
				calculate_number_of_goptima_per_acc_run_get_result(res, cpop, noacc, acclevels, pFunc);
				// Create data file.
				// /////////////////////////////////////////////////////////////////////////////////////////
				//
				std::fstream file;
				std::string filename = p.alldata;
				if (p.header && header){
					file.open(filename.c_str(), std::fstream::out);
					file << "#algorithm,problem,run,fes,time,";
					for (int iacc = 0; iacc < noacc; ++iacc) {
						file << "numGOacc" << iacc << ",";
					}
					file << "numberGO,popsize" << endl;
					header = false;
				} else {
					file.open(filename.c_str(), std::fstream::out | std::fstream::app);
				}
				int number_of_goptima = pFunc->get_no_goptima();
				file << p.algname << "," << "prob" << index << "," << irun << ",";
				file << data[ii][data[ii].size()-3] << "," << data[ii][data[ii].size()-2] << ",";
				for (int iii=0; iii < res.size(); iii++) {
					file << res[iii] << ",";
				}
				file << number_of_goptima <<"," << cpop.size() << endl;
				file.close();
			}
		}
		cout << endl;
		/* Clean up */
		delete pFunc;
	}

	return EXIT_SUCCESS;
}

int read_data_line(const std::vector<double> &data, TwoVectors &result, const int &dim) 
{
	vector<double> indi;
	vector<double> cost;
	int archive_action;
	int icounter = 0;
	for (icounter=0; icounter<dim; ++icounter) {
		indi.push_back(data[icounter]);
	}
	cost.push_back(data[icounter]); // fitness
	cost.push_back(data[++icounter]); // fes
	cost.push_back(data[++icounter]); // time
	archive_action = data[++icounter]; // action

	result = std::make_pair(indi, cost);
	return archive_action;
}


void calculate_number_of_goptima_per_acc_run_get_result(std::vector<double> &results, const Array2D &pop, 
		const int &noacc, const double *acclevels, CEC2013 *pFunc)
{
	int size = pop.size();
	/* for each accuracy level */
	for (int iacc = 0; iacc < noacc; ++iacc) {
		/* Calculate how many global optima are in the population */
		double accuracy = acclevels[iacc];
		std::vector< std::vector<double> > seeds; 
		if (size == 0) {
			results.push_back(0);
		} else {
			int count = how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho());
			results.push_back(count);
		}
	}
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

void load_file_with_archive_action(const std::string filename, 
		std::vector< std::vector<double> > &data, 
		const int &dim)
{
	std::ifstream file;
	file.open(filename.c_str(), std::fstream::in);
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		vector<double> all;
		double tmpdbl, fes, time;
		int archive_action;
		for (int i=0; i<dim; ++i) {
			iss >> tmpdbl;
			all.push_back(tmpdbl);
		}
		char tmpchar;
		iss >> tmpchar;		// = 
		iss >> tmpdbl;			// fitness 
		iss >> tmpchar;		// @
		iss >> fes;			// fes
		iss >> time;		// time
		iss >> archive_action;		// archive action

		all.push_back(tmpdbl);
		all.push_back(fes);
		all.push_back(time);
		all.push_back(archive_action);
		data.push_back(all);
	}
	// Close file
	file.close();
}

std::string operator+ (const string& a, const string& b) {
	std::string result(a);
	result += b;
	return result;
}
std::string vector_to_string(const std::vector<double> & vec)
{
	std::string result = "";
	for (int i=0; i< vec.size(); i++) {
		result += std::to_string(vec[i]);
	}
	return result;
}

