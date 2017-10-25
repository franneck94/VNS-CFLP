#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <numeric>
#include <memory>
#include <random>
#include <chrono>
#include <direct.h>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Helper.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                                HELPER FUNCTIONS                               */
/*********************************************************************************/

// Get Uniform Randomly Double Value
double get_random_double(const double lower, const double upper)
{
	random_device random_device;
	mt19937 engine{ random_device() };
	uniform_real_distribution<double> dist(lower, upper);

	return dist(engine);
}

// Get Uniform Randomly Integer Value
int get_random_int(const int lower, const int upper)
{
	random_device random_device;
	mt19937 engine{ random_device() };
	uniform_int_distribution<int> dist(lower, upper);

	return dist(engine);
}

// Prints Result of y_i and x_ij
void printFinalSolution(
	const vector<bool> &yi,
	const double &fx,
	const unsigned int &kMax,
	const unsigned int &tMax,
	const unsigned int &shaking_mode,
	const unsigned int &local_search_mode,
	const unsigned int &update_xij_mode,
	const unsigned int &init_mode,
	const string &directory,
	chrono::steady_clock::time_point &bvns_start,
	chrono::steady_clock::time_point &bvns_end)
{
	//                         (0=F, 1=B)      (0=G, 1=R)
	// kMax, kMin, tM,   shake,   local,   xij,  init,     
	string solution =
			to_string(kMax) + "_"
		+	to_string(tMax) + "_"
		+	to_string(shaking_mode);

	string full_directory = directory.substr(0, directory.find("/Testdaten/")) + "/Ergebnisse/"
		+ directory.substr(directory.find_last_of("/") + 1) + "/";
	cout << endl << "Soltuion: " << solution << "  f(x) = " + to_string(fx) << endl;

	//_mkdir(full_directory.c_str());

	ofstream ofile(full_directory + solution + ".txt", ofstream::out);

	ofile << "BEGIN PARAMS" << endl;
	ofile << "kMax: " + to_string(kMax) << endl;
	ofile << "tMax: " + to_string(tMax) << endl;

	if(shaking_mode == 0)
		ofile << "Shaking: " << "shakingKOperations" << endl;
	else if (shaking_mode == 1)
		ofile << "Shaking: " << "shakingKMaxOperations" << endl;
	else if (shaking_mode == 2)
		ofile << "Shaking: " << "shakingAssignments" << endl;
	else
		ofile << "Shaking: " << "shakingCosts" << endl;

	if(local_search_mode == 0)
		ofile << "LocalSearch: " << "First Improvment" << endl;
	else
		ofile << "LocalSearch: " << "Best Improvment" << endl;

	ofile << "Update Xij: " << "Vogel Approx." << endl;

	if (init_mode == 1)
		ofile << "Init: " << "RVNS" << endl;
	else
		ofile << "Init: " << "Random" << endl;

	ofile << "Soltuion: f(x) = " + to_string(fx) << endl;
	ofile << "Overall Time: " << chrono::duration_cast<chrono::microseconds>(bvns_end - bvns_start).count() << "s" << endl;
	ofile << "END PARAMS" << endl << endl;
	
	for (size_t i = 0; i != yi.size(); ++i)
	{
		if (yi[i] == 1)
		{
			ofile << "y[" << i << "] = 1" << endl;
		}
	}

	ofile.close();
}