#include "params.h"
using namespace std;

parameters::parameters() : 
	noruns(50), problems(20), alldata("alldata.csv"), header(false)
{ 
}

parameters::parameters(int &argc, char **argv) :
	noruns(50), problems(20), alldata("alldata.csv"), header(false)
{
	set_comand_line_values(argc, argv);
}


parameters::parameters(const parameters &tmp) :
	noruns(tmp.noruns),
	problems(tmp.problems),
	outfilename(tmp.outfilename),
	respath(tmp.respath),
	algname(tmp.algname),
	alldata(tmp.alldata),
	header(tmp.header)
{
}

parameters & parameters::operator=(const parameters &rhs)
{
	if (&rhs == this) {
		return *this;
	}
	noruns		= rhs.noruns;
	problems	= rhs.problems;
	outfilename	= rhs.outfilename;
	respath		= rhs.respath;
	algname		= rhs.algname;
	alldata		= rhs.alldata;
	header		= rhs.header;

	return *this;
}

void parameters::set_comand_line_values(int &argc, char **argv)
{
	char o;
	while ( (o = getopt(argc, argv, "r:p:o:P:a:D:H")) != -1 ) {
		switch (o) {
			case 'r':{ noruns = atoi(optarg); break; }
			case 'p':{ problems = atoi(optarg); break; }
			case 'o':{ outfilename = (optarg); break; }
			case 'P':{ respath = (optarg); break; }
			case 'a':{ algname = (optarg); break; }
			case 'D':{ alldata = (optarg); break; }
			case 'H':{ header = true; break; }
			default :{ print(); break; }
		}
	}
}

void parameters::print()
{
	cout << "Params..." << endl;
}
