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
#include "VNS.hpp"
#include "Solution.hpp"
#include "Helper.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*							INITIAL SOLUTIONS FOR VNS                            */
/*********************************************************************************/

// Initial Feasible Solution
void VNS::initialGreedy(
	vector <bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &fi,
	const vector<double> &dj,
	ublas::matrix<double> &cij,
	const double &sum_dj)
{
	unsigned int r_index = 0;

	// Closes Facilitys till Demand can be fullfilled
	while (canUpdateXij(bi, incumbent_solution, sum_dj))
	{
		r_index = select_randomly(incumbent_solution);
		incumbent_solution[r_index] = 0;
	}
	incumbent_solution[r_index] = 1;

	if (!updateXij(m_cij, m_dj, m_bi, incumbent_solution, m_update_xij_mode, m_best_fx))
	{
		m_best_fx = DBL_MAX;
		fill(incumbent_solution.begin(), incumbent_solution.end(), 1);
	}

	m_fx = m_best_fx;
}

// Inital Solution from RVNS
void VNS::initialRVNS(
	vector <bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &dj,
	ublas::matrix<double> &cij)
{
	unsigned int  k = 0, r_index = 0;

	// Closes Facilitys till Demand can be fullfilled
	while (canUpdateXij(bi, incumbent_solution, m_sum_dj))
	{
		r_index = select_randomly(incumbent_solution);
		incumbent_solution[r_index] = 0;
	}
	incumbent_solution[r_index] = 1;

	if (!updateXij(m_cij, m_dj, m_bi, incumbent_solution, m_update_xij_mode, m_fx))
	{
		m_fx = DBL_MAX;
		fill(incumbent_solution.begin(), incumbent_solution.end(), 1);
	}

	if (m_rMax > 0)
	{
		for (size_t i = 0; i != m_rMax; ++i)
		{
			// For all N_k Structures
			m_k = m_kMin - 1;

			while (m_k < m_kMax)
			{
				// Shake incumbent Solution
				m_perturbed_solution = shaking(incumbent_solution, bi, dj, cij, m_k, m_fx, m_shaking_mode);
				// Move or not
				neighborhoodChange(incumbent_solution, m_perturbed_solution, m_k, m_fx, m_best_fx, m_kMin);
			}
		}

		m_k = m_kMin - 1;
		m_fx = m_best_fx;
	}
}