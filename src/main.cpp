#include "runner.h"
#include "argparser.h"
#include "gen.h"


int main(int argc, char ** argv) {
	ArgParser parser(argc, argv);
	double beta = parser["beta"].as<double>();
	int n_atoms = 108, n_configs = parser["number"].as<int>();
	Generator g(n_atoms, parser);
	// targetEnergy(g1, n_configs);
	// tempSweep(g, n_configs, true);
	// cout << "\nFinding energy range...\n";
	// double emin = eMin(g, n_configs);
	// double emax = eMax(g, n_configs);
	// cout << "Min Coulomb energy: " << emin << endl;
	// cout << "Max Coulomb energy: " << emax << endl;
	// generateLattices(g, n_configs, emin, emax, 10, beta);
	// targetEnergy(g, n_configs, 540000, beta);
	equilTest(g, n_configs);
	return 0;
}