#include "runner.h"
#include "argparser.h"
#include "gen.h"
#include "stdlib.h"


void allocateRAMDisk() {
	system("sudo mkdir -p memory/");
	system("sudo mount -t tmpfs -o size=2048M tmpfs memory/");
}

void deallocateRAMDisk() {
	system("sudo umount memory/");
	system("sudo rm -r memory/");
}

int main(int argc, char ** argv) {
	allocateRAMDisk();
	ArgParser parser(argc, argv);
	double beta = parser["beta"].as<double>();
	double targetE = parser["etarget"].as<double>();
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
	// targetEnergy(g, n_configs, targetE, beta);
	equilTest(g, n_configs);
	deallocateRAMDisk();
	return 0;
}