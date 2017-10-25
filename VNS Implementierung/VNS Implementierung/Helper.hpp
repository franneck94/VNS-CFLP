#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <random>
#include <iterator>
#include <numeric>
#include <chrono>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                                HELPER FUNCTIONS                               */
/*********************************************************************************/

// Make Boost::Matrix from given Sizes m,n with Values v
template <typename T, typename F = ublas::row_major>
ublas::matrix<T, F> makeMatrix(size_t m, size_t n, vector<T> & v)
{
	ublas::unbounded_array<T> storage(m*n);
	copy(v.begin(), v.end(), storage.begin());
	return ublas::matrix<T>(m, n, storage);
}

// Print Vector Elements inline
template<typename T>
void printVector(const T& t)
{
	cout << "Vector: ";
	for (size_t i = 0; i != t.size(); ++i)
	{
		if (i == 0)
			cout << "( " << t[i] << " ,";
		else if (i == t.size() - 1)
			cout << " " << t[i] << " )";
		else
			cout << " " << t[i] << " ,";
	}
	cout << endl;
}

// Print Matrix Elements
template<typename T>
void printMatrix(const ublas::matrix<T> &m)
{
	cout << "Matrix: " << endl;
	for (size_t i = 0; i != m.size1(); ++i)
	{
		for (size_t j = 0; j != m.size2(); ++j)
		{
			if (j == 0)
				cout << "(" << m(i, j) << " , \t\t";
			else if(j == m.size2()-1)
				cout << " " << m(i, j) << ")";
			else
				cout << " " << m(i, j) << ", \t\t";
		}
		cout << endl;
	}
}

// Select Random index from Vector
template<typename T>
int select_randomly(const vector<T> &vec)
{
	if (vec.size() > 1)
	{
		random_device random_device;
		mt19937 engine{ random_device() };
		size_t upper = vec.size() - 1;
		uniform_int_distribution<int> dist(0, upper);

		return dist(engine);
	}

	return 0;
}

// Get Uniform Randomly Double or Integer Value
double get_random_double(const double lower, const double upper);
int get_random_int(const int lower, const int upper);

void printFinalSolution(
	const vector<bool> &yi,
	const double &fx,
	const unsigned int &kMax,
	const unsigned int &tMax,
	const unsigned int &shaking_mode,
	const unsigned int &local_search_mode,
	const unsigned int &update_xij_mode,
	const unsigned int &init_mode,
	const string &directory,
	chrono::steady_clock::time_point &bvns_start,
	chrono::steady_clock::time_point &bvns_end); // Print Solution Functions