#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "TestSet.hpp"
#include "Helper.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

/*********************************************************************************/
/*                                CLASS FOR TESTSET                              */
/*********************************************************************************/

TestSet::TestSet(const string &directory): m_directory(directory)
{
	// Read in Data from Files
	ReadInData(directory, m_dj, m_bi, m_fi, m_cij);

	// I, J Ranges
	m_locationNumber = m_fi.size();
	m_customerNumber = m_dj.size();
}

void TestSet::ReadInData(
	const string &directory,
	vector<double> &dj,
	vector<double> &bi,
	vector<double> &fi,
	ublas::matrix<double> &cij)
{
	string line = "";
	vector<string> myFiles = { "fi.txt", "bi.txt", "dj.txt" };
	vector <double> cache;
	double vector_value = 0.0;
	string::size_type sz;

	// Read in: d_j, f_i, b_i
	for (string file : myFiles)
	{
		ifstream infile(directory + '/' + file);

		for (string line; getline(infile, line); )
		{
			if (file == "dj.txt")
			{
				if (line.find_last_of(' ') != string::npos)
					vector_value = stod(line.substr(line.find_last_of(' ')).c_str());
				else if (!line.empty())
					vector_value = stod(line, &sz);
				dj.push_back(vector_value);
			}
			else if (file == "fi.txt")
			{
				if (line.find_last_of(' ') != string::npos)
					vector_value = stod(line.substr(line.find_last_of(' ')).c_str());
				else if (!line.empty())
					vector_value = stod(line, &sz);
				fi.push_back(vector_value);
			}
			else if (file == "bi.txt")
			{
				if (line.find_last_of(' ') != string::npos)
					vector_value = stod(line.substr(line.find_last_of(' ')).c_str());
				else if (!line.empty())
					vector_value = stod(line, &sz);
				bi.push_back(vector_value);
			}
		}

		infile.close();
	}

	// Read in C_ij
	ifstream infile(directory + '/' + "cij.txt");
	for (string line; getline(infile, line); )
	{
		if (line.find_last_of(' ') != string::npos)
			vector_value = stod(line.substr(line.find_last_of(' ')).c_str());
		else if (!line.empty())
			vector_value = stod(line, &sz);
		cache.push_back(vector_value);
	}

	cij = makeMatrix(bi.size(), dj.size(), cache);
	cache.clear();
	infile.close();
}