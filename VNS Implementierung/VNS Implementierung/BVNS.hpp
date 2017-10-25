#pragma once

#include "VNS.hpp"

using namespace std;

/*********************************************************************************/
/*                        BVNS CLASS, Inherit from VNS                          */
/*********************************************************************************/

class BVNS : public VNS
{
public:

	/*********************************************************************************/
	/*                                 CONSTRUCTOR                                   */
	/*********************************************************************************/

	BVNS(const string &directory, 
		const unsigned int &kMin,
		const unsigned int &kMax,
		const unsigned int &tMax,
		const unsigned int &shaking_mode,
		const unsigned int &local_Search_mode,
		const unsigned int &update_xij_mode,
		const unsigned int &initial_mode,
		const unsigned int &time); // Init BVNS Heuristic

	/*********************************************************************************/
	/*                FIRST OR BEST IMPROVMENT LOCAL SEARCH                          */
	/*********************************************************************************/

	vector<bool> localSearch(
		const vector<bool> &perturbed_solution,
		const vector<vector<bool>> &neighborhood_k,
		const unsigned int &k,
		double &fx,
		const unsigned int &local_search_mode); // Local Search Procedure

	vector<bool> localSearchFirstImprovment(
		const vector<bool> &perturbed_solution,
		const vector<vector<bool>> &neighborhood_k,
		const unsigned int &k,
		double &fx); // Search First Improvment in N_k(S)

	vector<bool> localSearchBestImprovment(
		const vector<bool> &perturbed_solution,
		const vector<vector<bool>> &neighborhood_k,
		const unsigned int &k,
		double &fx); // Search Best Improvment in N_k(S)

/*********************************************************************************/
/*                              MEMBER VARIABLES                                */
/*********************************************************************************/

private:
	vector<bool> m_local_solution; // Local Minimum Solution
};