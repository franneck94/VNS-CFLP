#include <iostream>
#include <string>
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

#include "VNS.hpp"
#include "Solution.hpp"
#include "Helper.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                   CONSTRUCTOR AND MEMBER FUNCTIONS OF VNS                     */
/*********************************************************************************/

VNS::VNS(const string &directory,
	const unsigned int &kMin,
	const unsigned int &kMax,
	const unsigned int &tMax,
	const unsigned int &shaking_mode,
	const unsigned int &local_search_mode,
	const unsigned int &update_xij_mode,
	const unsigned int &initial_mode,
	const unsigned int &time):
	// Init VNS and TestSet
		TestSet(directory), m_k(kMin - 1), m_kMax(kMax), m_kMin(kMin), m_tMax(tMax),
		m_shaking_mode(shaking_mode), m_local_search_mode(local_search_mode), m_update_xij_mode(update_xij_mode),
		m_incumbent_solution(m_locationNumber, 1), m_perturbed_solution(m_locationNumber, 1),
		m_fx(DBL_MAX), m_best_fx(DBL_MAX), m_sum_dj(accumulate(m_dj.begin(), m_dj.end(), 0.0)), 
		m_initial_mode(initial_mode), m_maxNH(pow(m_kMax, 2)),m_timeMax(time), m_rMax(0), 
		m_neighborhood_k(m_maxNH, vector<bool>(m_locationNumber, 1)), m_flow_tpl(0), m_initial_fx(DBL_MAX)
{
}

void VNS::initialSolution(
	vector <bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &fi,
	const vector<double> &dj,
	ublas::matrix<double> &cij,
	const double &sum_dj,
	const unsigned int &init_mode)
{
	if (init_mode == 0)
		initialGreedy(incumbent_solution, bi, fi, dj, cij, sum_dj);
	else
		initialRVNS(incumbent_solution, bi,  dj, cij);

	if (m_best_fx <= m_fx)
		m_initial_fx = m_best_fx;
	else
		m_initial_fx = m_fx;
}

// Move or not
void VNS::neighborhoodChange(
	vector<bool> &incumbent_solution,
	const vector<bool> &local_solution,
	unsigned int &k,
	const double &fx,
	double &best_fx,
	const double &kMin)
{
	if (fx < best_fx)
	{
		incumbent_solution = local_solution;
		k = kMin - 1;
		best_fx = fx;
	}
	else
		k++;
}

// Creates Perturbed Solution S'
vector<bool> VNS::shaking(
	const vector<bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &dj,
	ublas::matrix<double> &cij,
	const double &k,
	double &fx,
	const unsigned int &shaking_mode)
{
	if (shaking_mode == 0)
		return shakingKOperations(incumbent_solution, bi, m_dj, m_cij, m_k, fx);
	else if (shaking_mode == 1)
		return shakingKMaxOperations(incumbent_solution, bi, m_dj, m_cij, m_k, fx);
	else if (shaking_mode == 2)
		return shakingAssignments(incumbent_solution, bi, m_dj, m_cij, m_k, fx);
	else
		return shakingCosts(incumbent_solution, bi, m_dj, m_cij, m_k, fx);
}

// Getter f(x)
const double VNS::getFx() const
{
	return m_best_fx;
}

// Getter initial f(x)
const double VNS::getInitialFx() const
{
	return m_initial_fx;
}