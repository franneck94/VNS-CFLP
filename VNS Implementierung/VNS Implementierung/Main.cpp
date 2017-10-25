#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

#include "BVNS.hpp"
#include "Helper.hpp"

using namespace std;

int main(int argc, char **argv)
{
	unsigned int kMin = 1, kMax = 0, tMax = 0, shaking_mode = 0, 
		update_neighborhood_mode = 0, local_search_mode = 0, update_xij_mode = 0, 
		stopping_critera = 0, init_mode = 0, rMax = 0, time = 0, rMax_mode = 0;

	vector<double> vals(10, 0.0);
	vector<double> vals_init(10, 0.0);

	string test_set_directory = "C:/Users/Jan/Dropbox/Bachelorarbeit/Programm/Testdaten/cap";
	BVNS *bvns = nullptr;

	vector<string> small_sets = {"10", "11","12", "13", "14", "15"};
	vector<string> medium_sets = { "16", "17", "18", "19", "20", "21" };
	vector<string> large_sets = { "a", "b", "c", "d", "e", "f" };
	vector<string> larger_sets = {"22","23","24","25","26","27","28","29",
		"30","31","32","33","34","35","36", "37", "38", "39", "40", "41", "42", "43", 
		"44", "45", "46", "47", "48", "49", "50", "51"};
	vector<string> all_sets = { "10", "11","12", "13", "14", "15", "16", "17", "18", "19", 
		"20", "21","a", "b", "c", "d", "e", "f", "22","23","24","25","26","27","28","29",
		"30","31","32","33","34","35","36", "37", "38", "39", "40", "41", "42", "43",
		"44", "45", "46", "47", "48", "49", "50", "51" };

	//vector<string> todo_sets = { "16", "a", "22", "28", "36", "44" };
	vector<string> todo_sets = { "22" };

	for (auto set : todo_sets)
	{
		cout << "Open Testset: cap" << set << endl;

		local_search_mode = 1;
		shaking_mode = 0;
		update_xij_mode = 1;
		init_mode = 1;
		time = 3000;
		tMax = 40;
		kMax = 3;

		bvns = new BVNS(test_set_directory + set, kMin, kMax, tMax,
			shaking_mode, local_search_mode, update_xij_mode, init_mode, time);

		cout << endl << "init = " << bvns->getInitialFx() << endl;
	}

	cout << endl << endl << "Finished Computation!" << endl;
	getchar();
	return 0;
}

// Init: RVNS, Local: Best, Xij: Vogel, Time Cap: 30s
/*init_mode = 1;
local_search_mode = 1;
update_xij_mode = 1;
time = 30;
tMax = 40;
kMax = 3;*/