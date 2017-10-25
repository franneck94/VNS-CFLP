#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <tuple>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Solution.hpp"
#include "Helper.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

// Check if Solution is feasible
bool checkSolution(
	const vector<bool> &yi,
	const vector<tuple<int,int,double>> &flow,
	const unsigned int &flowNumber,
	const ublas::matrix<double> &cij,
	const vector<double> &bi,
	vector<double> dj,
	const double &sum_dj)
{
	double sum = 0.0;
	// Conditions are: forall j in J: sum of xij = 1
	//				   if xij !0 0 then yi = 1
	//				   sum_{j in J} xij * dj <= bi*yi, forall i in I
	for (size_t f = 0; f < flowNumber; f++)
	{
		if (yi[get<0>(flow[f])] == 0)
			return false;

		sum += get<2>(flow[f]);
		dj[get<1>(flow[f])] -= get<2>(flow[f]);
	}

	if (sum != sum_dj)
		return false;

	for (size_t j = 0; j != cij.size2(); ++j)
	{
		if (dj[j] != 0.0)
			return false;
	}
	
	return true;
}

// Obejtive Function of FLP
double f(
	const vector<bool> &yi,
	const vector<double> &fi,
	const double &transportation_cost)
{
	double z = 0.0;

	for (size_t i = 0; i != yi.size(); ++i)
		z += yi[i] * fi[i];

	z += transportation_cost;

	return z;
}

// Checks if Capacity can fullfill the Demand
bool canUpdateXij(
	const vector<double> &bi,
	const vector<bool> &yi,
	const double &sum_dj)
{
	if (inner_product(yi.begin(), yi.end(), bi.begin(), 0.0) < sum_dj)
		return false;
	else
		return true;
}