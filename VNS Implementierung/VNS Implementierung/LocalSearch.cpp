#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "BVNS.hpp"
#include "Helper.hpp"
#include "Solution.hpp"

using namespace std;

/*********************************************************************************/
/*                FIRST OR BEST IMPROVMENT LOCAL SEARCH                          */
/*********************************************************************************/

// Search First Improvment in N_k(S)
vector<bool> BVNS::localSearchFirstImprovment(
	const vector<bool> &perturbed_solution,
	const vector<vector<bool>> &neighborhood_k,
	const unsigned int &k,
	double &fx)
{
	int index = -1;
	unsigned int gap_mode = 0;
	double incumbent_value = fx;

	// Find first Imporvment in N_k(S')
	for (size_t i = 0; i != neighborhood_k.size(); ++i)
	{
		if (canUpdateXij(m_bi, neighborhood_k[i], m_sum_dj))
		{
			if (updateXij(m_cij, m_dj, m_bi, neighborhood_k[i], m_update_xij_mode, fx))
			{
				if (fx < incumbent_value)
				{
					index = i;
					break;
				}
			}
		}
	}
	// If Improvment was found
	if (index != -1)
		return neighborhood_k[index];
	// No Improvment found
	else
	{
		fx = incumbent_value;
		return perturbed_solution;
	}
}

// Search Best Improvment in N_k(S)
vector<bool> BVNS::localSearchBestImprovment(
	const vector<bool> &perturbed_solution,
	const vector<vector<bool>> &neighborhood_k,
	const unsigned int &k,
	double &fx)
{
	int index = -1;
	unsigned int gap_mode = 0;
	double incumbent_value = fx;
	double save_incumbent = fx;

	// Gets the best Solution in N_k(x) Neighborhood
	for (size_t i = 0; i != neighborhood_k.size(); ++i)
	{
		if (canUpdateXij(m_bi, neighborhood_k[i], m_sum_dj))
		{
			if (updateXij(m_cij, m_dj, m_bi, neighborhood_k[i], m_update_xij_mode, fx))
			{
				if (fx < incumbent_value)
				{
					index = i;
					incumbent_value = fx;
				}
			}
		}
	}

	// If Improvment was found
	if (index != -1)
	{
		fx = incumbent_value;
		return neighborhood_k[index];
	}
	// No Improvment found
	else
	{
		fx = save_incumbent;
		return perturbed_solution;
	}
}