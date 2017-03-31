#include "SimuInterface.h"
#include <algorithm>
using namespace std;

int main(int argc, char* argv[]){// arguments should be input file path
	cout << "Number of input files: " << argc-1 << endl;
	for (int i = 1; i < argc; ++i){
		cout << "Processing input file No." << i << " out of " << argc-1 << "..." << endl;
		SimuInterface simulator;
		bool success = false;

		string filename = string(argv[i]);
		if(filename.substr(filename.find_last_of(".") + 1) == "h5") 
		{
			if(filename.find("restart") != string::npos){
				success = simulator.import_restart_HDF5(argv[i]); 
			}
			else{
				success = simulator.import_HDF5(argv[i]);
			}
		} 
		else
		{
			cout << "Unrecogized input filename extension. " << endl;
		}
		
		if (success){ // return true if import is successful
			simulator.simulate();
			cout << "Input file No." << i << " out of " << argc-1 << " processed." << endl;
		}
	}
	cout << "The planet earth is blue and there's nothing I can do." << endl;
	return 0; 
};
