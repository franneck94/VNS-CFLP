#include <vector>

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "MODI.hpp"

namespace ublas = boost::numeric::ublas;
using namespace std;

namespace t_simplex
{
	/* DECLARATION OF GLOBALS */
	double ** _tsC = NULL;				// Cost matrix
	double _tsMaxW;						// Maximum of all weights

										// MODI Method
	double transportSimplex(
		const unsigned &update_xij_mode,
		ublas::matrix<double> &custom_cij,
		TsFlow *flowTable,
		unsigned int *flowSize,
		TsSignature *facility,
		TsSignature *customer, 
		const double &sum_dj,
		const double &sum_bi,
		const vector <double > &dj)
	{
		int m = facility->n, n = customer->n; // Matrix Sizes
		int i = 0, j = 0;
		double totalCost = 0.0, diff = 0.0;
		int *P1 = NULL, *P2 = NULL;

		TsBasic *basics = NULL;					///Array of basic variables. 
		bool **isBasic = NULL;						//Flag matrix. isBasic[i][j] is true there is flow between source i and sink j
		TsBasic **facBasics = NULL;					//Array of pointers to the first basic variable in each row
		TsBasic **cusBasics = NULL;					//Array of pointers to the first basic variable in each column
		double *demand = facility->weights;						//Array of sink demands
		double *capacity = NULL;					//Array of source supplies
		bool flag = false;

		// Equalize source and sink weights.
		diff = sum_bi - sum_dj;
		if (diff > 0.0) n++;
		_tsMaxW = sum_bi > sum_dj ? sum_bi : sum_dj;

		basics = new TsBasic[m + n];
		isBasic = new bool*[m];
		for (i = 0; i < m; ++i)
			isBasic[i] = NULL;
		for (i = 0; i < m; ++i)
		{
			isBasic[i] = new bool[n];
			for (j = 0; j < n; ++j)
				isBasic[i][j] = 0;
		}
		facBasics = new TsBasic*[m];
		for (i = 0; i < m; ++i)
			facBasics[i] = NULL;
		cusBasics = new TsBasic*[n];
		for (i = 0; i < n; ++i)
			cusBasics[i] = NULL;

		// Compute the cost matrix
		_tsC = new double*[m];
		for (i = 0; i < m; ++i)
			_tsC[i] = NULL;

		for (i = 0, P1 = facility->features; i < m; ++i, P1++) 
		{
			_tsC[i] = new double[n];					
			for (j = 0, P2 = customer->features; j < n; ++j, P2++)
			{
				if (i == facility->n || j == customer->n) 
					_tsC[i][j] = 0;			
				else
					_tsC[i][j] = custom_cij(i, j);
			}
		}

		capacity = new double[m];					//init the source array
		for (i = 0; i < facility->n; ++i) 
			capacity[i] = facility->weights[i];

		demand = new double[n];					//init the sink array
		for (i = 0; i < customer->n; ++i) 
			demand[i] = customer->weights[i];

		if (n != customer->n)
			demand[customer->n] = diff;

		// Find the initail basic feasible solution. Use either _initRussel or _initVogel
		_initVogel(capacity, demand, basics, facBasics, cusBasics, isBasic, m, n);

		// Enter the main pivot loop
		totalCost = _pivot(basics, facBasics, cusBasics, isBasic, m, n, dj);

		// Fill the Flow data structure
		TsBasic * basicPtr = NULL;
		TsFlow * flowPtr = flowTable;
		unsigned int counter = 0;

		if (flowTable != NULL) {
			for (i = 0; i < m + n; ++i) {
				basicPtr = basics + i;

				if (basicPtr != NULL)
				{
					if (basicPtr->i >= 0 && basicPtr->i < facility->n
						&& basicPtr->j >= 0 && basicPtr->j < customer->n
						&& isBasic[basicPtr->i][basicPtr->j]
						&& basicPtr->val != 0.0
						&& counter < m + n - 1)
					{
						flowPtr->to = basicPtr->j;
						flowPtr->from = basicPtr->i;
						flowPtr->amount = basicPtr->val;
						flowPtr++;
						counter++;
					}
				}
			}
		}
		if (flowSize != NULL)
		{
			*flowSize = (int)(flowPtr - flowTable);
		}

		for (i = 0; i < m; i++)
		{
			delete[] isBasic[i];
			isBasic[i] = NULL;
		}
		delete[] isBasic;
		isBasic = NULL;

		for (i = 0; i < m; i++)
		{
			delete[] _tsC[i];
			_tsC[i] = NULL;
		}
		delete[] _tsC;
		_tsC = NULL;

		delete[] facBasics;
		facBasics = NULL;
		delete[] cusBasics;
		cusBasics = NULL;
		delete[] basics;
		basics = NULL;

		return totalCost;
	}

	/*
	Main pivot loop.
	Pivots until the system is optimal and return the optimal transportation cost.
	*/
	double _pivot(
		TsBasic * basics,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int &m,
		int &n,
		const vector <double > &dj)
	{
		double * ui = NULL;
		double * vj = NULL;
		TsStone * stonePath = NULL;

		TsStone * spitra = NULL, *spitrb = NULL, *leaving = NULL;
		TsBasic * XP = NULL;
		TsBasic * basicsEnd = basics + m + n;
		TsBasic * entering = basicsEnd - 1;
		TsBasic dummyBasic;
		dummyBasic.i = -1;
		dummyBasic.j = 0;

		unsigned int i = 0, j = 0, lowI = 0, lowJ = 0, numPivots = 0;
		double objectiveValue = TSINFINITY, oldObjectiveValue = 0.0;
		double lowVal = 0.0;

		ui = new double[m];
		vj = new double[n];
		stonePath = new TsStone[m + n];

		while (1)
		{
			oldObjectiveValue = objectiveValue;
			objectiveValue = 0.0;

			for (XP = basics; XP != basicsEnd; XP++)
			{
				if (XP != entering)
				{
					//cout << endl << "i = " << XP->i << " j = " << XP->j << endl;
					if (XP->i >= 0 && XP->j >= 0 && XP->i < m && XP->j < dj.size()
						&& _tsC[XP->i][XP->j] != 0.0 && dj[XP->j] != 0.0)
					{
						objectiveValue += _tsC[XP->i][XP->j] * (XP->val / dj[XP->j]);
					}
				}
			}

			// Compute ui, vj
			stonePath[0].node = &dummyBasic;
			stonePath[0].prev = NULL;
			spitrb = _BFS(stonePath, facBasics, cusBasics, true);

			spitra = stonePath;
			vj[spitra->node->j] = 0;
			for (spitra++; spitra != spitrb; spitra++) {
				if (spitra->node->i == spitra->prev->node->i) {
					//node is in same row as parent
					vj[spitra->node->j] = _tsC[spitra->node->i][spitra->node->j] - ui[spitra->node->i];
				}
				else if (spitra->node->j == spitra->prev->node->j) {
					ui[spitra->node->i] = _tsC[spitra->node->i][spitra->node->j] - vj[spitra->node->j];
				}
			}

			// find Theta
			lowVal = 0.0;
			for (i = 0; i < m; ++i)
				for (j = 0; j < n; ++j)
					if (!isBasic[i][j] && _tsC[i][j] - ui[i] - vj[j] < lowVal)
					{
						lowVal = _tsC[i][j] - ui[i] - vj[j];
						lowI = i;
						lowJ = j;
					}

			if (lowVal >=  0.0 || (oldObjectiveValue - objectiveValue) < TSPIVOTLIMIT)
			{
				delete[] ui;
				delete[] vj;
				delete[] stonePath;
				//std::cout << "Pivots: " << numPivots << "\t";
				return objectiveValue;
			}

			// Add the entering variable to stone path
			entering->i = lowI;
			entering->j = lowJ;
			isBasic[lowI][lowJ] = 1;
			entering->val = 0;
			entering->nextFac = facBasics[lowI];
			if (facBasics[lowI] != NULL) facBasics[lowI]->prevFac = entering;
			entering->nextCus = cusBasics[lowJ];
			if (cusBasics[lowJ] != NULL) cusBasics[lowJ]->prevCus = entering;
			facBasics[lowI] = entering;
			entering->prevFac = facBasics[lowI];
			cusBasics[lowJ] = entering;
			entering->prevCus = cusBasics[lowJ];
			stonePath[0].node = entering;
			stonePath[0].prev = NULL;

			// Use breadth-first search to find a loop of basics.
			spitra = spitrb = _BFS(stonePath, facBasics, cusBasics);
			lowVal = TSINFINITY;
			bool add = false;

			// Find the lowest flow along the loop (leaving variable)
			do
			{
				if (!add && spitrb->node->val < lowVal)
				{
					leaving = spitrb;
					lowVal = spitrb->node->val;
				}
				add = !add;
			} while (spitrb = spitrb->prev);

			add = false;
			spitrb = spitra;

			// Alternately increase and decrease flow along the loop
			do
			{
				if (add)
					spitrb->node->val += lowVal;
				else
					spitrb->node->val -= lowVal;
				add = !add;
			} while (spitrb = spitrb->prev);

			i = leaving->node->i;
			j = leaving->node->j;
			isBasic[i][j] = 0;

			if (facBasics[i] == leaving->node)
			{
				facBasics[i] = leaving->node->nextFac;
				facBasics[i]->prevFac = NULL;
			}
			else
			{
				leaving->node->prevFac->nextFac = leaving->node->nextFac;
				if (leaving->node->nextFac != NULL)
					leaving->node->nextFac->prevFac = leaving->node->prevFac;
			}

			if (cusBasics[j] == leaving->node)
			{
				cusBasics[j] = leaving->node->nextCus;
				cusBasics[j]->prevCus = NULL;
			}
			else
			{
				leaving->node->prevCus->nextCus = leaving->node->nextCus;
				if (leaving->node->nextCus != NULL)
					leaving->node->nextCus->prevCus = leaving->node->prevCus;
			}

			entering = leaving->node;
			numPivots++;
		}
	}

	TsStone * _BFS(
		TsStone * stoneTree,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool complete)
	{
		bool column = true;
		int jumpoffset = 0;
		TsBasic * bitr;
		TsStone * sitra = &stoneTree[0], *sitrb = &stoneTree[1];
		do {
			if (column)
			{
				for (bitr = cusBasics[sitra->node->j]; bitr != NULL; bitr = bitr->nextCus)
				{
					if (bitr != sitra->node) 
					{
						sitrb->node = bitr;
						sitrb->prev = sitra;
						sitrb++;
					}
				}
			}
			else 
			{
				for (bitr = facBasics[sitra->node->i]; bitr != NULL; bitr = bitr->nextFac)
				{
					if (bitr != sitra->node)
					{
						sitrb->node = bitr;
						sitrb->prev = sitra;
						sitrb++;
					}
				}
			}

			sitra++;
			if (sitra == sitrb) //no cycle found and no cycles in tree
				return sitra;

			if (sitra->node->i == sitra->prev->node->i)
				column = true;
			else
				column = false;

			// cycle found
			if (!complete && sitra->node->i == stoneTree[0].node->i
				&& sitra->node->j != stoneTree[0].node->j  && column == false)
				return sitra;
		} while (1);
	}

	/**********************
	Vogel's initialization method
	**********************/
	void _initVogel(
		double *S,
		double *D,
		TsBasic * basicsEnd,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int m,
		int n)
	{
		int i, j;
		TsVogPen *srcPens = NULL;
		TsVogPen *snkPens = NULL;
		TsVogPen *pitra, *pitrb;  //iterators
		TsVogPen *maxPen;
		TsVogPen srcPenHead, snkPenHead;
		bool maxIsSrc;
		double lowVal;

		srcPens = new TsVogPen[m];
		snkPens = new TsVogPen[n];

		srcPenHead.next = pitra = srcPens;
		for (i = 0; i < m; i++)
		{
			pitra->i = i;
			pitra->next = pitra + 1;
			pitra->prev = pitra - 1;
			pitra->one = pitra->two = 0;
			pitra->oneCost = pitra->twoCost = TSINFINITY;
			pitra++;
		}
		(--pitra)->next = NULL;
		srcPens[0].prev = &srcPenHead;

		snkPenHead.next = pitra = snkPens;
		for (i = 0; i < n; i++) 
		{
			pitra->i = i;
			pitra->next = pitra + 1;
			pitra->prev = pitra - 1;
			pitra->one = pitra->two = 0;
			pitra->oneCost = pitra->twoCost = TSINFINITY;
			pitra++;
		}
		(--pitra)->next = NULL;
		snkPens[0].prev = &snkPenHead;


		for (pitra = srcPenHead.next, i = 0; pitra != NULL; pitra = pitra->next, i++)
		{
			for (pitrb = snkPenHead.next, j = 0; pitrb != NULL; pitrb = pitrb->next, j++)
			{
				//initialize Source Penalties;
				addPenalty(pitra, _tsC[i][j], j);
				addPenalty(pitrb, _tsC[i][j], i);
			}
		}

		while (srcPenHead.next != NULL && snkPenHead.next != NULL)
		{
			maxIsSrc = true;
			for (maxPen = pitra = srcPenHead.next; pitra != NULL; pitra = pitra->next)
				if ((pitra->twoCost - pitra->oneCost) > (maxPen->twoCost - maxPen->oneCost))
					maxPen = pitra;

			for (pitra = snkPenHead.next; pitra != NULL; pitra = pitra->next)
				if ((pitra->twoCost - pitra->oneCost) > (maxPen->twoCost - maxPen->oneCost))
				{
					maxPen = pitra;
					maxIsSrc = false;
				}

			if (maxIsSrc)
			{
				i = maxPen->i;
				j = maxPen->one;
			}
			else 
			{
				j = maxPen->i;
				i = maxPen->one;
			}

			if (D[j] - S[i] > _tsMaxW * TSEPSILON || (srcPenHead.next->next != NULL && fabs(S[i] - D[j]) < _tsMaxW * TSEPSILON))
			{
				//delete source
				lowVal = S[i];
				maxPen = srcPens + i;
				maxPen->prev->next = maxPen->next;
				if (maxPen->next != NULL)
					maxPen->next->prev = maxPen->prev;

				for (pitra = snkPenHead.next; pitra != NULL; pitra = pitra->next) {
					if (pitra->one == i || pitra->two == i) {
						pitra->oneCost = TSINFINITY;
						pitra->twoCost = TSINFINITY;
						for (pitrb = srcPenHead.next; pitrb != NULL; pitrb = pitrb->next)
							addPenalty(pitra, _tsC[pitrb->i][pitra->i], pitrb->i);
					}
				}
			}
			else 
			{
				//delete sink
				lowVal = D[j];
				maxPen = snkPens + j;
				maxPen->prev->next = maxPen->next;
				if (maxPen->next != NULL)
					maxPen->next->prev = maxPen->prev;

				for (pitra = srcPenHead.next; pitra != NULL; pitra = pitra->next)
				{
					if (pitra->one == j || pitra->two == j) {
						pitra->oneCost = TSINFINITY;
						pitra->twoCost = TSINFINITY;
						for (pitrb = snkPenHead.next; pitrb != NULL; pitrb = pitrb->next)
							addPenalty(pitra, _tsC[pitra->i][pitrb->i], pitrb->i);
					}
				}
			}

			S[i] -= lowVal;
			D[j] -= lowVal;

			isBasic[i][j] = 1;
			basicsEnd->val = lowVal;
			basicsEnd->i = i;
			basicsEnd->j = j;

			basicsEnd->nextCus = cusBasics[j];
			if (cusBasics[j] != NULL) cusBasics[j]->prevCus = basicsEnd;
			basicsEnd->nextFac = facBasics[i];
			if (facBasics[i] != NULL) facBasics[i]->prevFac = basicsEnd;

			facBasics[i] = basicsEnd;
			basicsEnd->prevCus = NULL;
			cusBasics[j] = basicsEnd;
			basicsEnd->prevFac = NULL;

			basicsEnd++;

		}
		delete[] srcPens;
		delete[] snkPens;
	}
}