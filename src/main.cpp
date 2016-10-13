#include <iostream>
#include <thread>
#include <boost/program_options.hpp> 

#include "gen.h"
#include "argparser.h"

using namespace std;

namespace po = boost::program_options;

int main(int argc, char ** argv) {
	ArgParser parser(argc, argv);
	int n_atoms = 108;
	Generator g(n_atoms, parser);
	cout << g.totalCoulombEnergy() << "\n";
	return 0;
}