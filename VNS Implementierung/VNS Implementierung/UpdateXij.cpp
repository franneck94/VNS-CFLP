#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <random>
#include <iterator>
#include <numeric>

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "VNS.hpp"
#include "MODI.hpp"
#include "Solution.hpp"
#include "Helper.hpp"

using namespace std;
using namespace t_simplex;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                            UPDATE XIJ PROCEDURES                              */
/*********************************************************************************/

// CFLP with Modi procedure (gets optimal solution in deterministic time)
bool VNS::updateXij(
	ublas::matrix<double> &cij,
	const vector<double> &dj,
	const vector<double> &bi,
	const vector<bool> &yi,
	const unsigned int &update_xij_mode,
	double &fx)
{
	if (!canUpdateXij(bi, yi, m_sum_dj))
		return false;

	unsigned int open_facilities = 0, i_out = 0, flow_vars = 0;
	double transportation_cost = 0.0, gap = 0.0, sum_bi = 0.0;
	bool check = false;

	open_facilities = accumulate(yi.begin(), yi.end(), 0);

	ublas::matrix<double> custom_cij(open_facilities, m_customerNumber, 0.0);

	int *customers = new int[m_customerNumber];
	double *demand = new double[m_customerNumber];
	int *facilities = new int[open_facilities];
	double *capacity = new double[open_facilities];

	for (size_t j = 0; j != cij.size2(); ++j)
	{
		customers[j] = j;
		demand[j] = dj[j];
	}

	for (size_t i = 0; i != cij.size1(); ++i)
	{
		if (yi[i] == 1)
		{
			facilities[i_out] = i;
			capacity[i_out] = bi[i];
			sum_bi += bi[i];

			ublas::matrix_row<ublas::matrix<double>> row_custom_cij(custom_cij, i_out);
			ublas::matrix_row<ublas::matrix<double>> row_cij(cij, i);
			row_custom_cij = row_cij;

			i_out++;
		}
	}

	// Signature of Facilities and Customers
	TsSignature *facility = new TsSignature(open_facilities, facilities, capacity);
	TsSignature *customer = new TsSignature(m_customerNumber, customers, demand);
	// Save Stepping Stone Path
	TsFlow *flow = new TsFlow[open_facilities + m_customerNumber - 1];

	// Result value
	transportation_cost = t_simplex::transportSimplex(update_xij_mode, custom_cij, flow, 
		&flow_vars, facility, customer ,m_sum_dj, sum_bi, dj);

	// f(x) and Return Value
	fx = f(yi, m_fi, transportation_cost);

	// Flow Xij for checkSolution
	m_flow_tpl.resize(flow_vars);
	for (size_t i = 0; i < flow_vars; ++i)
		m_flow_tpl[i] = make_tuple(facilities[flow[i].from], flow[i].to, flow[i].amount);

	check = checkSolution(yi, m_flow_tpl, flow_vars, cij, bi, dj, m_sum_dj);

	delete[] capacity;
	capacity = NULL;
	delete[] demand;
	demand = NULL;
	delete[] facilities;
	facilities = NULL;
	delete[] customers;
	customers = NULL;
	delete facility;
	facility = NULL;
	delete customer;
	customer = NULL;

	return check;
}