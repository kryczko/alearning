#ifndef __ARGPARSER__
#define __ARGPARSER__

#include <boost/program_options.hpp> 
#include <string> 

namespace po = boost::program_options;

using namespace std;

class ArgParser {
    private:
        po::variables_map vm;
        po::options_description desc;

    public:

        po::variable_value operator[](const string option) {
            if (this->vm.count(option)) {
                return this->vm[option];
            } else {
                cout << "ERROR: " << option << " was not found.\n";
                exit(-1);
            }
        }

        ArgParser(int argc, char ** argv) {
            this->desc.add_options()
            ("help", "Generate FCC Al lattices of user defined sizes. Example usage: $ ./generator --defects 1 -dopants Ti 5 O 3 --number 100")
            ("defects", po::value<int>()->default_value(0), "Number of defects to introduce into the lattice.")
            ("dopants", po::value<vector<string>>()->multitoken()->default_value(vector<string>(), ""), "Dopants type followed by number of dopant atoms to insert into lattice.")
            ("number", po::value<int>()->default_value(10000), "Number of configurations to generate.")
            ("etarget", po::value<double>()->default_value(1000000.), "Target energy when sampling lattices.")
            ("seed", po::value<int>()->default_value(13), "Random seed.")
            ("beta", po::value<double>()->default_value(1.0), "Inverse temperature used in the Bolzmann factor.");
            po::store(po::parse_command_line(argc, argv, this->desc), this->vm);
            po::notify(this->vm);

            if (this->vm.count("help")) {
                cout << this->desc << "\n";
                exit(-1);
            }
        }
};

#endif