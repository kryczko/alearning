#include "runner.h"
#include "argparser.h"
#include "gen.h"


int main(int argc, char ** argv) {
	ArgParser parser(argc, argv);
	int n_atoms = 108, n_configs = parser["number"].as<int>();
	Generator g(n_atoms, parser);
	// targetEnergy(g1, n_configs);
	// tempSweep(g1, n_configs, true);
	equilTest(g, n_configs);
	return 0;
}