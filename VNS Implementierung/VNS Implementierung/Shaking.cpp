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
/*                                SHAKING PROCEDURES                             */
/*********************************************************************************/

// Creates Perturbed Solution S'
vector<bool> VNS::shakingKOperations(
	const vector<bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &dj,
	ublas::matrix<double> &cij,
	const double &k,
	double &fx)
{
	size_t last = 0, lastS = 0, counter = 0, i1 = 0, i2 = 0;
	vector<bool> perturbed_solution = incumbent_solution;
	double rnd = 0.0;

	vector<int> I_minus_S(incumbent_solution.size());
	vector<int> S(incumbent_solution.size());

	// I \ S := all closed Facilities
	last = 0;
	lastS = 0;
	for (size_t index = 0; index != perturbed_solution.size(); ++index)
		if (perturbed_solution[index] == false)
			I_minus_S[last++] = index;
		else
			S[lastS++] = index;
	I_minus_S.erase(I_minus_S.begin() + last, I_minus_S.end());
	S.erase(S.begin() + lastS, S.end());

	do
	{
		// Uniformed Random Value
		rnd = get_random_double(0.0, 1.0);

		// i1 element S
		if (S.size() != 0)
			i1 = S[select_randomly(S)];
		else
			i1 = 0;

		// i2 element I \ S
		if (I_minus_S.size() != 0)
			i2 = I_minus_S[select_randomly(I_minus_S)];
		else
			i2 = 0;

		// Drop
		if (rnd <= 0.2)
		{
			perturbed_solution[i1] = 0;
			I_minus_S.push_back(i1);
		}
		// Add
		else if (rnd >= 0.8)
		{
			perturbed_solution[i2] = 1;
			I_minus_S.erase(remove(
				I_minus_S.begin(), I_minus_S.end(), i2), I_minus_S.end());
		}
		// Swap
		else
		{
			perturbed_solution[i1] = 0;
			perturbed_solution[i2] = 1;

			I_minus_S.push_back(i1);
			I_minus_S.erase(remove(
				I_minus_S.begin(), I_minus_S.end(), i2), I_minus_S.end());
		}

		counter++;
	} while (counter < k); // Repeat k-Times

	// If its in Shaking-Phase and not in creating N_k(S) 
	if (fx != -1.0)
	{
		if (!canUpdateXij(bi, perturbed_solution, m_sum_dj))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		// Check Solution if its not feasible return incumbent
		if (!updateXij(cij, dj, bi, perturbed_solution, m_update_xij_mode, fx))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		return perturbed_solution;
	}

	return perturbed_solution;
}

// Creates Perturbed Solution S'
vector<bool> VNS::shakingKMaxOperations(
	const vector<bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &dj,
	ublas::matrix<double> &cij,
	const double &k,
	double &fx)
{
	size_t last = 0, lastS = 0, counter = 0, i1 = 0, i2 = 0;
	vector<bool> perturbed_solution = incumbent_solution;
	double rnd = 0.0;

	vector<int> I_minus_S(incumbent_solution.size());
	vector<int> S(incumbent_solution.size());

	// I \ S := all closed Facilities
	last = 0;
	lastS = 0;
	for (size_t index = 0; index < perturbed_solution.size(); ++index)
		if (perturbed_solution[index] == false)
			I_minus_S[last++] = index;
		else
			S[lastS++] = index;
	I_minus_S.erase(I_minus_S.begin() + last, I_minus_S.end());
	S.erase(S.begin() + lastS, S.end());

	do
	{
		// Uniformed Random Value
		rnd = get_random_double(0.0, 1.0);

		// i1 element S
		if (S.size() != 0)
			i1 = S[select_randomly(S)];
		else
			i1 = 0;

		// i2 element I \ S
		if (I_minus_S.size() != 0)
			i2 = I_minus_S[select_randomly(I_minus_S)];
		else
			i2 = 0;

		// Drop
		if (rnd <= 0.2)
		{
			perturbed_solution[i1] = 0;
			I_minus_S.push_back(i1);
		}
		// Add
		else if (rnd >= 0.8)
		{
			perturbed_solution[i2] = 1;
			I_minus_S.erase(remove(
				I_minus_S.begin(), I_minus_S.end(), i2), I_minus_S.end());
		}
		// Swap
		else
		{
			perturbed_solution[i1] = 0;
			perturbed_solution[i2] = 1;

			I_minus_S.push_back(i1);
			I_minus_S.erase(remove(
				I_minus_S.begin(), I_minus_S.end(), i2), I_minus_S.end());
		}		

		counter++;
	} while (counter < m_kMax); // Repeat k-Times

	// If its in Shaking-Phase and not in creating N_k(S) 
	if (fx != -1.0)
	{
		if (!canUpdateXij(bi, perturbed_solution, m_sum_dj))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		// Check Solution if its not feasible return incumbent
		if (!updateXij(cij, dj, bi, perturbed_solution, m_update_xij_mode, fx))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		return perturbed_solution;
	}

	return perturbed_solution;
}

// Shaking Method from Kratica et al.
vector<bool> VNS::shakingAssignments(
	const vector<bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &dj,
	ublas::matrix<double> &cij,
	const double &k,
	double &fx)
{
	unsigned int nh_mode = get_random_int(0, 2), i_out = 0, rand_indx = 0;
	double yi_sum = 0.0;
	size_t last = 0, lastS = 0, counter = 0, i1 = 0, i2 = 0;
	vector<bool> perturbed_solution = incumbent_solution;
	double rnd = 0.0;

	vector<int> I_minus_S(incumbent_solution.size());
	vector<int> S(incumbent_solution.size());

	// I \ S := all closed Facilities
	last = 0;
	lastS = 0;
	for (size_t index = 0; index < perturbed_solution.size(); ++index)
		if (perturbed_solution[index] == false)
			I_minus_S[last++] = index;
		else
			S[lastS++] = index;
	I_minus_S.erase(I_minus_S.begin() + last, I_minus_S.end());
	S.erase(S.begin() + lastS, S.end());
	
	// Swap
	if (nh_mode == 2)
	{
		for (size_t c = 0; c != k; ++c)
		{
			i1 = S[select_randomly(S)];
			i2 = I_minus_S[select_randomly(I_minus_S)];

			perturbed_solution[i1] = 0;
			perturbed_solution[i2] = 1;
		}
	}
	// Close k-Min/Max and open Random
	else
	{
		vector<double> assignments(m_locationNumber, 0);

		// Save all assignments
		for (size_t i = 0; i != m_flow_tpl.size(); ++i)
		{
			if (perturbed_solution[get<0>(m_flow_tpl[i])] == 1)
			{
				assignments[get<0>(m_flow_tpl[i])] += get<2>(m_flow_tpl[i]);
			}
		}

		// Index Sort Assignment list
		yi_sum = accumulate(incumbent_solution.begin(), incumbent_solution.end(), 0);
		vector<unsigned int> sorted_index;

		for (size_t indx = 0; indx != perturbed_solution.size(); ++indx)
		{
			if (perturbed_solution[indx] == 1)
			{
				sorted_index.push_back(indx);
			}
		}

		// sorted_index Index list
		sort(begin(sorted_index), end(sorted_index),
			[&](int i1, int i2) { return assignments[i1] < assignments[i2]; });

		if (nh_mode == 1)
		{
			for (size_t c = 0; c != k; ++c)
			{
				if (sorted_index.size() != 0 && c < perturbed_solution.size())
				{
					// Close k-Max
					if (c < sorted_index.size())
						perturbed_solution[sorted_index[c]] = 0;
					// Close k-Min
					if (sorted_index.size() - 1 - c >= 0)
						perturbed_solution[sorted_index[c]] = 0;
				}
			}
			for (size_t c = 0; c != k; ++c)
			{
				rand_indx = select_randomly(perturbed_solution);
				perturbed_solution[rand_indx] = 1;
			}
		}
		else
		{
			for (size_t c = 0; c != k; ++c)
			{
				if (sorted_index.size() != 0 && c < perturbed_solution.size())
				{
					// Close k-Min
					if (sorted_index.size() - 1 - c >= 0)
						perturbed_solution[sorted_index[c]] = 0;
				}
			}
			for (size_t c = 0; c != k; ++c)
			{
				rand_indx = select_randomly(perturbed_solution);
				perturbed_solution[rand_indx] = 1;
			}
		}	
	}

	// If its in Shaking-Phase and not in creating N_k(S) 
	if (fx != -1.0)
	{
		if (!canUpdateXij(bi, perturbed_solution, m_sum_dj))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		// Check Solution if its not feasible return incumbent
		if (!updateXij(cij, dj, bi, perturbed_solution, m_update_xij_mode, fx))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		return perturbed_solution;
	}

	return perturbed_solution;
}

// Modified Shaking Method from Kratica et al.
vector<bool> VNS::shakingCosts(
	const vector<bool> &incumbent_solution,
	const vector<double> &bi,
	const vector<double> &dj,
	ublas::matrix<double> &cij,
	const double &k,
	double &fx)
{
	unsigned int nh_mode = get_random_int(0, 2), i_out = 0, rand_indx = 0;
	double yi_sum = 0.0;
	size_t last = 0, lastS = 0, counter = 0, i1 = 0, i2 = 0;
	vector<bool> perturbed_solution = incumbent_solution;
	double rnd = 0.0;

	vector<int> I_minus_S(incumbent_solution.size());
	vector<int> S(incumbent_solution.size());

	// I \ S := all closed Facilities
	last = 0;
	lastS = 0;
	for (size_t index = 0; index < perturbed_solution.size(); ++index)
		if (perturbed_solution[index] == false)
			I_minus_S[last++] = index;
		else
			S[lastS++] = index;
	I_minus_S.erase(I_minus_S.begin() + last, I_minus_S.end());
	S.erase(S.begin() + lastS, S.end());

	// Normal Operations
	if (nh_mode == 1)
	{
		for (size_t c = 0; c != k; ++c)
		{
			i1 = S[select_randomly(S)];
			i2 = I_minus_S[select_randomly(I_minus_S)];

			perturbed_solution[i1] = 0;
			perturbed_solution[i2] = 1;
		}
	}
	// Close k-Min/Max and open Random
	else
	{
		vector<double> costs(m_locationNumber, 0);

		// Save all Costs per Facility
		for (size_t i = 0; i != cij.size1(); ++i)
		{
			ublas::matrix_row<ublas::matrix<double>> row_cij(cij, i);

			// For every xij link
			for (size_t f = 0; f != m_flow_tpl.size(); ++f)
			{
				if (get<0>(m_flow_tpl[f]) == i)
					costs[i] += row_cij[get<1>(m_flow_tpl[f])] * get<2>(m_flow_tpl[f]);
			}
		}

		// Index Sort Assignment list
		yi_sum = accumulate(perturbed_solution.begin(), perturbed_solution.end(), 0);
		vector<unsigned int> sorted_index;

		for (size_t indx = 0; indx != perturbed_solution.size(); ++indx)
		{
			if (perturbed_solution[indx] == 1)
			{
				sorted_index.push_back(indx);
			}
		}

		// sorted_index Index list
		sort(begin(sorted_index), end(sorted_index),
			[&](int i1, int i2) { return costs[i1] < costs[i2]; });

		if (nh_mode == 1)
		{
			for (size_t c = 0; c != k; ++c)
			{
				if (sorted_index.size() != 0 && c < perturbed_solution.size())
				{
					// Close k-Max
					if (c < sorted_index.size())
						perturbed_solution[sorted_index[c]] = 0;
				}
			}
			for (size_t c = 0; c != k; ++c)
			{
				do
				{
					rand_indx = select_randomly(perturbed_solution);
				} while (perturbed_solution[rand_indx] == 0);

				perturbed_solution[rand_indx] = 1;
			}

			for (size_t c = 0; c != k; ++c)
			{
				if (sorted_index.size() != 0 && c < perturbed_solution.size())
				{
					// Close k-Min
					if (sorted_index.size() - 1 - c >= 0)
						perturbed_solution[sorted_index[c]] = 0;
				}
			}
			for (size_t c = 0; c != k; ++c)
			{
				do
				{
					rand_indx = select_randomly(perturbed_solution);
				} while (perturbed_solution[rand_indx] == 0);

				perturbed_solution[rand_indx] = 1;
			}
		}
		// nh-gap_mode == 2
		else
		{
			for (size_t c = 0; c != k; ++c)
			{
				if (sorted_index.size() != 0 && c < perturbed_solution.size())
				{
					// Close k-Max
					if (c < sorted_index.size())
						perturbed_solution[sorted_index[c]] = 0;
				}
			}
			for (size_t c = 0; c != k; ++c)
			{
				rand_indx = select_randomly(perturbed_solution);
				perturbed_solution[rand_indx] = 1;
			}
		}	
	}

	// If its in Shaking-Phase and not in creating N_k(S) 
	if (fx != -1.0)
	{
		if (!canUpdateXij(bi, perturbed_solution, m_sum_dj))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		// Check Solution if its not feasible return incumbent
		if (!updateXij(cij, dj, bi, perturbed_solution, m_update_xij_mode, fx))
		{
			fx = m_best_fx;
			return incumbent_solution;
		}

		return perturbed_solution;
	}

	return perturbed_solution;
}