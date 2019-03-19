/******************************************************************************
 * Version: 1.0.1
 * Last modified on: 19 March, 2019 
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
void load_file(const std::string, std::vector< std::vector<double> > &, std::vector< std::vector<double> > &, std::vector< std::vector<double> > &, const int &);
void load_file(const std::string, std::vector< std::vector<double> > &, const int &);
void load_file_with_archive_action(const std::string filename, std::vector< std::vector<double> > &data, const int &dim);
int read_data_line(const std::vector<double> &data, TwoVectors &result, const int &dim) ;
int read_line(const std::string &line, TwoVectors &result, const int &dim) ;
void load_file_keep_last_archive(const std::string filename, std::vector< std::vector<double> > &pop, std::vector< std::vector<double> > &costs, const int &dim);
void init_memory(Array3D &results, Array2D &PR, Array2D &SR, 
		const int &problems, const int &noruns, const int &noacc);
void gather_stats(Array3D &results, Array2D &PR, Array2D &SR,
		const int &problems, const int &noruns, const int &noacc);
void create_latex_table(Array2D &PR, Array2D &SR, const std::string &outfilename, 
		double *acclevels, const std::string &algname, const int &problems, const int &noacc);
void re_evaluate_solutions(std::vector<double> &fitnesses, const Array2D &pop, CEC2013 *pFunc);
void re_evaluate_solutions(std::vector<double> &fitnesses, const Array2D &pop, Array2D &data, CEC2013 *pFunc);
void calculate_distances_from_closest_optimum(Array2D &all_dists, const Array2D &pop, const Array2D &goptima);
void calculate_number_of_goptima_per_acc_run(Array3D &results, const Array2D &pop, 
		const int &index, const int &irun, const int &noacc, const double *acclevels, CEC2013 *pFunc);
void calculate_number_of_goptima_per_acc_run_get_result(std::vector<double> &results, const Array2D &pop, 
		const int &noacc, const double *acclevels, CEC2013 *pFunc);
void calculate_PR_SR_per_problem(Array3D &results, Array2D &PR, Array2D &SR,
		const int &problems, const int &noruns, const int &noacc, const int &index, CEC2013 *pFunc);
void create_data_files(Array2D &PR, Array2D &SR, const std::string &outfilename, 
		double *acclevels, const std::string &algname, const int &problems, const int &noacc);
void create_data_file(const parameters &p, const int &index, const int &irun, const Array2D &pop, const std::vector<double> &fitnesses, const Array2D &costs, const Array2D &goptima, CEC2013 *pFunc, const Array2D &all_dists, Array3D &results, const int &noacc, const double *acclevels);
void create_data_file_dynamic_res(const parameters &p, const int &index, const int &irun, const Array2D &pop, const std::vector<double> &fitnesses, const Array2D &costs, const Array2D &goptima, CEC2013 *pFunc, const Array2D &all_dists, Array2D &results, const int &noacc, const double *acclevels);

std::string int_to_string_with_leading_zeros(int no, size_t orders_of_magnitute)
{
	std::ostringstream ss; ss << setw(orders_of_magnitute) << setfill('0') << no; return ss.str();
}
struct VectorComparatorFEs {
	bool operator()(const std::vector<double>& a, const std::vector<double> & b) {
		return (a[a.size()-2] < b[b.size()-2] );
	}
};
struct VectorComparatorTime {
	bool operator()(const std::vector<double>& a, const std::vector<double> & b) {
		return (a[a.size()]-1 < b[b.size()-1] );
	}
};
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

			// incrementally add points
			Array2D res_all;
			// create archive
			archive myarchive;
			for (unsigned int ii=0; ii< data.size(); ii++) {
				myarchive.clear();
				for (unsigned int jj=0; jj<=ii; jj++) {
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
				res_all.push_back(res);
			}
			// here I need the last pop, costs and fitnesses
			load_file_keep_last_archive(tmpfilename, pop, costs, dim);
			// re-evaluate solutions
			re_evaluate_solutions(fitnesses, pop, data, pFunc);
			
			// calculate distances from the closest optimum (decision space)
			calculate_distances_from_closest_optimum(all_dists, pop, goptima);

			// Calculate how many global optima are in the population
//			calculate_number_of_goptima_per_acc_run(results, pop, index, irun, noacc, acclevels, pFunc);

			// Create data file.
			create_data_file_dynamic_res(p, index, irun, pop, fitnesses, costs, goptima, pFunc, all_dists, res_all, noacc, acclevels);
			//return 0;
		}
		cout << endl;
		// calculate PR, SR per problem
		// calculate_PR_SR_per_problem(results, PR, SR, problems, noruns, noacc, index, pFunc);
		/* Clean up */
		delete pFunc;
	}
//	gather_stats(results, PR, SR, problems, noruns, noacc);
//	create_data_files(PR, SR, outfilename, acclevels, algname, problems, noacc);

	/* Create the PR/SR latex table */
	//create_latex_table(PR, SR, outfilename, acclevels, algname, problems, noacc);

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


void load_file_keep_last_archive(const std::string filename, 
		std::vector< std::vector<double> > &pop, 
		std::vector< std::vector<double> > &costs, 
		const int &dim)
{
	archive myarchive;
	std::ifstream file;
	file.open(filename.c_str(), std::fstream::in);
	std::string line;
	int linecount = 1 ;
	while (std::getline(file, line)) {
		TwoVectors result;
		int action = read_line(line, result, dim);
		if (action == 0) {
			// empty the archive and 
			myarchive.clear();
			// add the new solution
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
				cerr << "WARNING: I could not remove solution at line " << linecount << " because it was not in the archive." << endl;
			}
		} else {
			cerr << "ERROR: " << filename << " contains incorrect archive action index: " << action << " at line: "<< linecount << endl;
			exit(-1);
		}
		linecount++;
	}
	// Close file
	file.close();
	// add info from archive into pop/costs
	pop.clear(); costs.clear();
	for (archive::iterator it = myarchive.begin(); it != myarchive.end(); ++it) {
		pop.push_back((it->second).first);
		costs.push_back((it->second).second);
	}
}

int read_line(const std::string &line, TwoVectors &result, const int &dim) 
{
	std::istringstream iss(line);
	vector<double> indi;
	vector<double> cost;
	int archive_action;
	double tmp, fes, time;
	for (int i=0; i<dim; ++i) {
		iss >> tmp;
		indi.push_back(tmp);
	}
	char tmpchar;
	iss >> tmpchar;		// = 
	iss >> tmp;			// fitness 
	iss >> tmpchar;		// @
	iss >> fes;			// fes
	iss >> time;		// time
	iss >> archive_action;		// archive action
	cost.push_back(fes);
	cost.push_back(time);
	result = std::make_pair(indi, cost);
	return archive_action;
}

void create_data_file_dynamic_res(const parameters &p, const int &index, const int &irun, const Array2D &pop, const std::vector<double> &fitnesses, const Array2D &costs, const Array2D &goptima, CEC2013 *pFunc, const Array2D &all_dists, Array2D &results, const int &noacc, const double *acclevels)
{
	std::fstream file;
	std::string filename = p.alldata;
	if (p.header && header){
		file.open(filename.c_str(), std::fstream::out);
		file << "#algorithm,problem,run,counter,";
		file << "fitness,fitnessGO,fes,time,dist,idxdist,";
		for (int iacc = 0; iacc < noacc; ++iacc) {
			file << "numGOacc" << iacc << ",";
		}
		file << "numberGO" << endl;
		header = false;
	} else {
		file.open(filename.c_str(), std::fstream::out | std::fstream::app);
	}
	double fitGO = pFunc->evaluate(goptima[0]);
	int number_of_goptima = pFunc->get_no_goptima();
	for (int i=0; i < (int)pop.size(); i++) {
		file << p.algname << "," << "prob" << index << "," << irun << "," << (i+1) << ",";
		file << fitnesses[i] << "," << fitGO << "," << costs[i][0] << "," << costs[i][1] << ",";
		int min_dist_idx = all_dists[i][goptima.size()];
		double mindist = all_dists[i][min_dist_idx];
	   	file << mindist << "," << min_dist_idx << ",";
		for (int iacc = 0; iacc < noacc; ++iacc) {
			file << results[i][iacc] << ",";
		}
		file << number_of_goptima << endl;
	}

	file.close();
}

void create_data_file(const parameters &p, const int &index, const int &irun, const Array2D &pop, const std::vector<double> &fitnesses, const Array2D &costs, const Array2D &goptima, CEC2013 *pFunc, const Array2D &all_dists, Array3D &results, const int &noacc, const double *acclevels)
{
	std::fstream file;
	std::string filename = p.alldata;
	if (p.header && header){
		file.open(filename.c_str(), std::fstream::out);
		file << "#algorithm,problem,run,counter,";
		file << "fitness,fitnessGO,fes,time,dist,idxdist,";
		for (int iacc = 0; iacc < noacc; ++iacc) {
			file << "numGOacc" << iacc << ",";
		}
		file << "numberGO" << endl;
		header = false;
	} else {
		file.open(filename.c_str(), std::fstream::out | std::fstream::app);
	}
	double fitGO = pFunc->evaluate(goptima[0]);
	int number_of_goptima = pFunc->get_no_goptima();
	for (int i=0; i < (int)pop.size(); i++) {
		file << p.algname << "," << "prob" << index << "," << irun << "," << (i+1) << ",";
		file << fitnesses[i] << "," << fitGO << "," << costs[i][0] << "," << costs[i][1] << ",";
		int min_dist_idx = all_dists[i][goptima.size()];
		double mindist = all_dists[i][min_dist_idx];
	   	file << mindist << "," << min_dist_idx << ",";
		for (int iacc = 0; iacc < noacc; ++iacc) {
			file << results[index-1][irun][iacc] << ",";
		}
		file << number_of_goptima << endl;
	}

	file.close();
}

void calculate_PR_SR_per_problem(Array3D &results, Array2D &PR, Array2D &SR,
		const int &problems, const int &noruns, const int &noacc, const int &index, CEC2013 *pFunc)
{
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
void calculate_number_of_goptima_per_acc_run(Array3D &results, const Array2D &pop, 
		const int &index, const int &irun, const int &noacc, const double *acclevels, CEC2013 *pFunc)
{
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

void calculate_distances_from_closest_optimum(Array2D &all_dists, const Array2D &pop, const Array2D &goptima)
{
	// calculate distances from the closest optimum (decision space)
	all_dists.clear();
	for (unsigned i=0; i < pop.size(); i++) {
		std::vector<double> dists;
		for (unsigned j=0; j < goptima.size(); j++) {
			double distance = 0.0;
			assert(pop[0].size() == goptima[0].size());
			for (unsigned k=0; k < goptima[j].size(); k++) {
				distance += (goptima[j][k] - pop[i][k]) * 
					(goptima[j][k] - pop[i][k]);
			}
			dists.push_back( sqrt(distance) );
		}
		int min_pos = std::distance(dists.begin(), std::min_element(dists.begin(),dists.end()) );
		dists.push_back( min_pos );
		all_dists.push_back(dists);
	}
}

void re_evaluate_solutions(std::vector<double> &fitnesses, const Array2D &pop, Array2D &data, CEC2013 *pFunc)
{
	fitnesses.clear();
	double fit;
	for (unsigned i=0; i < pop.size(); i++) {
		fit = pFunc->evaluate(pop[i]);
		fitnesses.push_back(fit);
		data[i][pFunc->get_dimension()] = fit;
	}
}

void re_evaluate_solutions(std::vector<double> &fitnesses, const Array2D &pop, CEC2013 *pFunc)
{
	fitnesses.clear();
	double fit;
	for (unsigned i=0; i < pop.size(); i++) {
		fit = pFunc->evaluate(pop[i]);
		fitnesses.push_back(fit);
	}
}

void create_data_files(Array2D &PR, Array2D &SR, const std::string &outfilename, 
		double *acclevels, const std::string &algname, const int &problems, const int &noacc)
{
	std::fstream filePR, fileSR;
	std::string filenamePR = outfilename + "_PR.csv";
	std::string filenameSR = outfilename + "_SR.csv";
	filePR.open(filenamePR.c_str(), std::fstream::out);
	fileSR.open(filenameSR.c_str(), std::fstream::out);

	/* for each problem */
	for (int i=0; i<problems; ++i) {
		/* for each accuracy level */
		for (int iacc = 0; iacc < noacc; ++iacc) {
			if (iacc == noacc-1) {
				filePR << PR[i][iacc];
				fileSR << SR[i][iacc];
			} else {
				filePR << PR[i][iacc] << ",";
				fileSR << SR[i][iacc] << ",";
			}
		}
		filePR << endl;
		fileSR << endl;
	}

	// Close file
	filePR.close();
	fileSR.close();
}

void create_latex_table(Array2D &PR, Array2D &SR, const std::string &outfilename, 
		double *acclevels, const std::string &algname, const int &problems, const int &noacc)
{
	const int rowblocks = 4;
	const int algs_per_rowblock = 5;
	std::fstream file;
	std::string tmpf = outfilename + ".tex";
	file.open(tmpf.c_str(), std::fstream::out);
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

void load_file(const std::string filename, 
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

		all.push_back(tmpdbl);
		all.push_back(fes);
		all.push_back(time);
		data.push_back(all);
	}
	// Close file
	file.close();
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

void load_file(const std::string filename, 
		std::vector< std::vector<double> > &pop, 
		std::vector< std::vector<double> > &costs, 
		std::vector< std::vector<double> > &data, 
		const int &dim)
{
	std::ifstream file;
	file.open(filename.c_str(), std::fstream::in);
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		vector<double> indi;
		vector<double> cost;
		vector<double> all;
		double tmp, fes, time;
		for (int i=0; i<dim; ++i) {
			iss >> tmp;
			indi.push_back(tmp);
			all.push_back(tmp);
		}
		pop.push_back(indi);
		char tmpchar;
		iss >> tmpchar;		// = 
		iss >> tmp;			// fitness 
		iss >> tmpchar;		// @
		iss >> fes;			// fes
		iss >> time;		// time
		cost.push_back(fes);
		cost.push_back(time);
		costs.push_back(cost);

		all.push_back(fes);
		all.push_back(time);
		data.push_back(all);
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
std::string vector_to_string(const std::vector<double> & vec)
{
	std::string result;
	for (int i=0; i< vec.size(); i++) {
		result += std::to_string(vec[i]);
	}
	return result;
}

