#ifndef __HELPER__
#define __HELPER__

#include <unordered_map>
#include <string>
#include <cmath>

typedef std::unordered_map<std::string, int> str_to_int_dict;

const double pbc = 8.568;
const double half_pbc = pbc / 2.0;
const str_to_int_dict coulomb_mapping (
	{
		{"H",1},
		{"He",1},
		{"Li",1},
		{"Be",1},
		{"B",1},
		{"C",1},
		{"N",1},
		{"O",1},
		{"F",1},
		{"Ne",1},
		{"Na",1},
		{"Mg",1},
		{"Al",1},
		{"Si",1},
		{"P",1},
		{"S",1},
		{"Cl",1},
		{"Ar",1},
		{"K",1},
		{"Ca",1},
		{"Sc",1},
		{"Ti",1},
		{"V",1},
		{"Cr",1},
		{"Mn",1},
		{"Fe",1},
		{"Co",1},
		{"Ni",1},
		{"Cu",1},
		{"Zn",1},
		{"Ga",1},
		{"Ge",1},
		{"As",1},
		{"Se",1},
		{"Br",1},
		{"Kr",1},
		{"Rb",1},
		{"Sr",1},
		{"Y",1},
		{"Zr",1},
		{"Nb",1},
		{"Mo",1}
	}
);


double pbc_wrap(double val) {
	if (val > half_pbc) {
		return val -= half_pbc;
	}
	return val;
}

double r(std::vector<double> a, std::vector<double> b) {
	double x = pbc_wrap(abs(a[0] - b[0]));
	double y = pbc_wrap(abs(a[1] - b[1]));
	double z = pbc_wrap(abs(a[2] - b[2]));
	return sqrt(x*x + y*y + z*z);
}

#endif