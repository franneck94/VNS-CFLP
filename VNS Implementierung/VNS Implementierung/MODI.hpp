#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#define TSINFINITY       1e20
#define TSEPSILON        1e-6
#define TSPIVOTLIMIT     0.00

using namespace std;
namespace ublas = boost::numeric::ublas;

namespace t_simplex {

	/* DECLARATION OF DATA TYPES */
	enum TsError { TsErrBadInput };

	//TsSignature is used for inputting the source and sink signatures

	class TsSignature
	{
	public:
		unsigned int n;									// Number of features in the signature 
		int *features;							// Pointer to the features vector 
		double *weights;						// Pointer to the weights of the features 
		TsSignature(unsigned int nin, int *fin, double * win) :n(nin), features(fin), weights(win) {};
	};

	//TsFlow is used for outputting the final flow table
	typedef struct TsFlow
	{
		int from;								// Feature number in signature 1 
		int to;									// Feature number in signature 2 
		double amount;							// Amount of flow from signature1.features[from] to signature2.features[to]
	} TsFlow;

	// TsBasic is used for 2D lists, allowing for easy navigation of the basic variables 
	typedef struct TsBasic
	{
		int i, j;
		double val;
		TsBasic *nextCus, *prevCus;				//next/previous node in the column
		TsBasic *nextFac, *prevFac;				//next/previous node in the row
	} TsBasic;

	// TsStone is used for _BFS
	typedef struct TsStone
	{
		struct TsStone *prev;
		struct TsBasic *node;
	} TsStone;

	// TsVogPen is used for 1D lists in _initVogel
	typedef struct TsVogPen
	{
		int i;
		struct TsVogPen *next, *prev;
		int one, two;
		double oneCost, twoCost;
	} TsVogPen;

	// Helper function for _initVogel
	inline void addPenalty(
		TsVogPen * pitr,
		double &cost,
		int &i)
	{
		if (pitr != NULL)
		{
			if (cost < pitr->oneCost)
			{
				pitr->twoCost = pitr->oneCost;
				pitr->two = pitr->one;
				pitr->oneCost = cost;
				pitr->one = i;
			}
			else if (cost < pitr->twoCost)
			{
				pitr->twoCost = cost;
				pitr->two = i;
			}
		}
	}

	/* DECLARATIONS */
	double transportSimplex(
		const unsigned &update_xij_mode,
		ublas::matrix<double> &custom_cij,
		TsFlow *flowTable,
		unsigned int *flowSize,
		TsSignature *facility,
		TsSignature *customer,
		const double &sum_dj,
		const double &sum_bi,
		const vector <double > &dj);

	double _pivot(
		TsBasic * basics,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int &m,
		int &n,
		const vector <double > &dj);

	TsStone * _BFS(
		TsStone * stoneTree,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool complete = false);

	void _initVogel(
		double *S,
		double *D,
		TsBasic * basicsEnd,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int m,
		int n);

	void _initNW(
		double *bi,
		double *dj,
		TsBasic * basicsEnd,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int m,
		int n);

	void _initLCM(
		double *bi,
		double *dj,
		TsBasic * basicsEnd,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int m,
		int n,
		ublas::matrix<unsigned int> &cij_column_sorted,
		int *facilities);
}
