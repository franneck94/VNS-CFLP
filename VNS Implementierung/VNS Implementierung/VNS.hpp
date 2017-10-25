#pragma once

#include "TestSet.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                      VNS CLASS, Inherit from TestSet                         */
/*********************************************************************************/

class VNS: public TestSet
{
public:

	/*********************************************************************************/
	/*                    CONSTRUCTOR AND NEIGHBORHOOD CHANGE                        */
	/*********************************************************************************/

	VNS(const string &directory,
		const unsigned int &kMin,
		const unsigned int &kMax,
		const unsigned int &tMax,
		const unsigned int &shaking_mode,
		const unsigned int &local_Search_mode,
		const unsigned int &update_xij_mode,
		const unsigned int &initial_mode,
		const unsigned int &time); // Can set all Options

	void neighborhoodChange(
		vector<bool> &incumbent_solution,
		const vector<bool> &local_solution,
		unsigned int &k,
		const double &fx,
		double &best_fx,
		const double &kMin); // Move or not

	/*********************************************************************************/
	/*                       UPDATE NEIGHBORHOOD PROCEDURES                          */
	/*********************************************************************************/

	void initialSolution(
		vector <bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &fi,
		const vector<double> &dj,
		ublas::matrix<double> &cij,
		const double &sum_dj,
		const unsigned int &init_mode); // Initial Feasible Solutio

	void initialGreedy(
		vector <bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &fi,
		const vector<double> &dj,
		ublas::matrix<double> &cij,
		const double &sum_dj); // Initial Feasible Solution from Random

	void initialRVNS(
		vector <bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &dj,
		ublas::matrix<double> &cij); // Initial Feasible Solution from RVNS

	/*********************************************************************************/
	/*                       UPDATE NEIGHBORHOOD PROCEDURE                           */
	/*********************************************************************************/

	void updateNeighborhoods(
		vector<vector<bool>> &neighborhood_k,
		const vector<bool> &perturbed_solution,
		const vector<double> &bi,
		const double &sum_dj,
		const unsigned int &shaking_mode,
		const double &act_k); // Updates Neighborhood from Solution S

	/*********************************************************************************/
	/*                                SHAKING PROCEDURES                             */
	/*********************************************************************************/

	vector<bool> shaking(
		const vector<bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &dj,
		ublas::matrix<double> &cij,
		const double &k,
		double &fx,
		const unsigned int &shaking_mode); // Creates Perturbed Solution S'

	vector<bool> shakingKOperations(
		const vector<bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &dj,
		ublas::matrix<double> &cij,
		const double &k,
		double &fx); // Uniformed Add/Remove/Swap Moves

	vector<bool> shakingKMaxOperations(
		const vector<bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &dj,
		ublas::matrix<double> &cij,
		const double &k,
		double &fx); // Uniformed Add/Remove/Swap Moves

	vector<bool> shakingAssignments(
		const vector<bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &dj,
		ublas::matrix<double> &cij,
		const double &k,
		double &fx); // Shaking Method from Kratica et al.

	vector<bool> shakingCosts(
		const vector<bool> &incumbent_solution,
		const vector<double> &bi,
		const vector<double> &dj,
		ublas::matrix<double> &cij,
		const double &k,
		double &fx); // Modified Shaking Method from Kratica et al.

	/*********************************************************************************/
	/*                            UPDATE XIJ PROCEDURES                              */
	/*********************************************************************************/

	bool updateXij(
		ublas::matrix<double> &cij,
		const vector<double> &dj,
		const vector<double> &bi,
		const vector<bool> &yi,
		const unsigned int &update_xij_mode,
		double &fx); // Update CFLP Matrix with MODI

	/*********************************************************************************/
	/*									   GETTER                                    */
	/*********************************************************************************/
	const double getFx() const;
	const double getInitialFx() const;

/*********************************************************************************/
/*                              MEMBER VARIABLES                                */
/*********************************************************************************/

protected:
	// VNS Parameters
	const unsigned int m_kMin;						// Neighborhood min. Num.
	const unsigned int m_kMax;						// Neighborhood max. Num.
	unsigned int m_k;								// Actual k-th Neighborhood
	const unsigned int m_maxNH;						// Max. Numb. of N_k(S)
	// VNS Setting Parameters
	const unsigned int m_tMax;						// Stopping Criteria e.g. max. Num. Iterations
	const unsigned int m_rMax;						// Max Iterations for RVNS init
	const unsigned int m_initial_mode;				// Inital Feasible Solution Mode
	const unsigned int m_shaking_mode;				// Shaking Mode		
	const unsigned int m_local_search_mode;			// First or Best Improvemnt
	const unsigned int m_update_xij_mode;			// Update Xij Mode (Minimum Searches)
	const unsigned m_timeMax;						// Stopping Criteria for Time against Gurobi
	// Soluton Parameters
	double m_initial_fx;							// Initial Solution Value
	double m_fx;									// Incumbent Solution Value
	double m_best_fx;								// Best Solution Value
	const double m_sum_dj;							// Sum of Demand
	vector<bool> m_incumbent_solution;				// S
	vector<bool> m_perturbed_solution;				// S'
	vector<vector<bool>> m_neighborhood_k;			// N_k(S)
	vector<tuple<int, int, double>> m_flow_tpl;		// Flow Tuple for xij Assignment
};