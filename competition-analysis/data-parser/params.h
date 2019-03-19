#ifndef parameters_h
#define parameters_h

#include <string>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <getopt.h>

class parameters
{ 
public:
	int noruns;
	int problems;
	std::string outfilename, respath, algname, alldata;
	bool header;

	parameters();
	parameters(int &argc, char **argv);
	parameters(const parameters &tmp);
	parameters& operator=(const parameters &rhs);
	void print();
private:
	void set_comand_line_values(int &argc, char **argv);
};
#endif
