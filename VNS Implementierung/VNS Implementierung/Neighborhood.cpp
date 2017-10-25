#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <numeric>
#include <memory>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "VNS.hpp"
#include "Helper.hpp"
#include "Solution.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                       UPDATE NEIGHBORHOOD PROCEDURES                          */
/*********************************************************************************/

// Creates or Updates Neighborhood Structure from Sol. S
void VNS::updateNeighborhoods(
	vector<vector<bool>> &neighborhood_k,
	const vector<bool> &perturbed_solution,
	const vector<double> &bi,
	const double &sum_dj,
	const unsigned int &shaking_mode,
	const double &k)
{
	double fx = -1;
	
	if (shaking_mode == 0)
	{
		// Iterate over all Neighborhood Solutions
		for (size_t n = 0; n < m_maxNH; ++n)
		{
			neighborhood_k[n] = shakingKOperations(perturbed_solution, bi, m_dj, m_cij, m_k, fx);
		}
	}
	else if(shaking_mode == 1)
	{
		// Iterate over all Neighborhood Solutions
		for (size_t n = 0; n < m_maxNH; ++n)
		{
			neighborhood_k[n] = shakingKMaxOperations(perturbed_solution, bi, m_dj, m_cij, m_k, fx);
		}
	}
	else if (shaking_mode == 2)
	{
		// Iterate over all Neighborhood Solutions
		for (size_t n = 0; n < m_maxNH; ++n)
		{
			neighborhood_k[n] = shakingAssignments(perturbed_solution, bi, m_dj, m_cij, m_k, fx);
		}
	}
	else
	{
		// Iterate over all Neighborhood Solutions
		for (size_t n = 0; n < m_maxNH; ++n)
		{
			neighborhood_k[n] = shakingCosts(perturbed_solution, bi, m_dj, m_cij, m_k, fx);
		}
	}
}