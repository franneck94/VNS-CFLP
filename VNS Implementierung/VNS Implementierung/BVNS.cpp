#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <numeric>
#include <memory>
#include <chrono>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "TestSet.hpp"
#include "BVNS.hpp"
#include "VNS.hpp"
#include "Solution.hpp"
#include "Helper.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                  CONSTRUCTOR AND MEMBER FUNCTIONS OF BVNS                     */
/*********************************************************************************/

// Init BVNS Heuristic
BVNS::BVNS(const string &directory, 
	const unsigned int &kMin,
	const unsigned int &kMax,
	const unsigned int &tMax,
	const unsigned int &shaking_mode,
	const unsigned int &local_Search_mode,
	const unsigned int &update_xij_mode,
	const unsigned int &initial_mode,
	const unsigned int &time):
		VNS(directory, kMin, kMax, tMax, shaking_mode, 
			local_Search_mode, update_xij_mode, initial_mode, time),
		m_local_solution(m_locationNumber, 1)
{
	chrono::steady_clock::time_point timeStart = chrono::high_resolution_clock::now();
	initialSolution(m_incumbent_solution, m_bi, m_fi, m_dj, m_cij, m_sum_dj, m_initial_mode);

	unsigned int counter = 0, t = 0;
	bool stop = false;
	
	for (size_t i = 0; i != tMax && !stop; ++i)
	{
		// For all N_k Structures
		m_k = m_kMin - 1;
		while (m_k < m_kMax && !stop)
		{
			// Shake incumbent Solution
			m_perturbed_solution = shaking(m_incumbent_solution,
				m_bi, m_dj, m_cij, m_k, m_fx, m_shaking_mode);
			// Update Neighborhood Structure N_k(S')
			updateNeighborhoods(m_neighborhood_k, m_perturbed_solution, m_bi, 
				m_sum_dj, m_shaking_mode, m_k);
			// Generate Local Solution
			m_local_solution = localSearch(m_perturbed_solution, 
				m_neighborhood_k, m_k, m_fx, m_local_search_mode);
			// Move or not
			neighborhoodChange(m_incumbent_solution, m_local_solution, m_k, m_fx, m_best_fx, m_kMin);

			// Check if maximum COmputation Time is exeeded
			chrono::steady_clock::time_point bvns_step = chrono::high_resolution_clock::now();
			double time_end = (double)chrono::duration_cast<chrono::microseconds>(bvns_step - timeStart).count();
			if (time_end >= m_timeMax * 1000000.0)
			{
				cout << "Time Limit:" << m_timeMax << " exceeded." << endl;
				stop = true;
			}
		}
	}
	
	chrono::steady_clock::time_point timeEnd = chrono::high_resolution_clock::now();

	// Compute Optimal MODI Value
	updateXij(m_cij, m_dj, m_bi, m_incumbent_solution, m_update_xij_mode, m_best_fx);

	// Save Solution and Timestamps to File
	printFinalSolution(m_incumbent_solution, m_best_fx, m_kMax, m_tMax, m_shaking_mode, 
		m_local_search_mode, m_update_xij_mode, m_initial_mode, m_directory, timeStart, timeEnd);
}

vector<bool> BVNS::localSearch(
	const vector<bool> &perturbed_solution,
	const vector<vector<bool>> &neighborhood_k,
	const unsigned int &k,
	double &fx,
	const unsigned int &local_search_mode)
{
	if (local_search_mode == 0)
		return localSearchFirstImprovment(perturbed_solution, neighborhood_k, k, fx);
	else
		return localSearchBestImprovment(perturbed_solution, neighborhood_k, k, fx);
}