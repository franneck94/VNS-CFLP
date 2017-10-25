#include <vector>
#include <algorithm>

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "MODI.hpp"
#include "Helper.hpp"

namespace ublas = boost::numeric::ublas;
using namespace std;

namespace t_simplex
{
	/**********************
	Vogel's initialization method
	**********************/
	void _initVogel(
		vector<double> &bi,
		vector<double> dj,
		TsBasic * basicsEnd,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int m,
		int n,
		ublas::matrix<double> &cij,
		double &_tsMaxW)
	{
		int i, j;
		TsVogPen *srcPens = NULL;
		TsVogPen *snkPens = NULL;
		TsVogPen *pitra = NULL, *pitrb = NULL;  //iterators
		TsVogPen *maxPen = NULL;
		TsVogPen srcPenHead, snkPenHead;
		bool maxIsSrc = false;
		double lowVal = 0.0;

		try
		{
			srcPens = new TsVogPen[m];
			snkPens = new TsVogPen[n];
		}
		catch (std::bad_alloc)
		{
			delete[] srcPens;
			delete[] snkPens;
			throw;
		}

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
				addPenalty(pitra, cij(i, j), j);
				addPenalty(pitrb, cij(i, j), i);
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

			if (dj[j] - bi[i] > _tsMaxW * TSEPSILON 
				|| (srcPenHead.next->next != NULL
				&& fabs(dj[i] - bi[j]) < _tsMaxW * TSEPSILON))
			{
				//delete source
				lowVal = bi[i];
				maxPen = srcPens + i;
				maxPen->prev->next = maxPen->next;
				if (maxPen->next != NULL)
					maxPen->next->prev = maxPen->prev;

				for (pitra = snkPenHead.next; pitra != NULL; pitra = pitra->next)
				{
					if (pitra->one == i || pitra->two == i)
					{
						pitra->oneCost = TSINFINITY;
						pitra->twoCost = TSINFINITY;
						for (pitrb = srcPenHead.next; pitrb != NULL; pitrb = pitrb->next)
							addPenalty(pitra, cij(pitrb->i, pitra->i), pitrb->i);
					}
				}
			}
			else
			{
				//delete sink
				lowVal = dj[j];
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
							addPenalty(pitra, cij(pitra->i, pitrb->i), pitrb->i);
					}
				}
			}

			bi[i] -= lowVal;
			dj[j] -= lowVal;

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
		srcPens = NULL;
		delete[] snkPens;
		snkPens = NULL;
	}

	/**********************
	NW-Corner initialization method
	**********************/
	void _initNW(
		double *bi,
		double *dj,
		TsBasic * basicsEnd,
		TsBasic ** facBasics,
		TsBasic ** cusBasics,
		bool ** isBasic,
		int m,
		int n)
	{
		unsigned int i = 0, j = 0;
		double lowVal = 0.0;

		while (i < m && j < n)
		{
			// More capacity than Demand
			if (bi[i] >= dj[j] && dj[j] != 0.0)
			{
				lowVal = dj[j];
				bi[i] -= dj[j];
				dj[j] = 0.0;

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

				// Test on Degeneration
				if (bi[i] == 0.0 && dj[j] == 0.0  && i < m - 1 && j < n - 1)
				{
					isBasic[i][j + 1] = 1;
					basicsEnd->val = lowVal;
					basicsEnd->i = i;
					basicsEnd->j = j + 1;

					basicsEnd->nextCus = cusBasics[j + 1];
					if (cusBasics[j + 1] != NULL) cusBasics[j + 1]->prevCus = basicsEnd;
					basicsEnd->nextFac = facBasics[i];
					if (facBasics[i] != NULL) facBasics[i]->prevFac = basicsEnd;

					facBasics[i] = basicsEnd;
					basicsEnd->prevCus = NULL;
					cusBasics[j + 1] = basicsEnd;
					basicsEnd->prevFac = NULL;

					basicsEnd++;

					i++;
					j++;
				}
				else
				{
					// Skip to next Customer
					j++;
				}
			}
			// Less Capacity than Demand
			else if (bi[i] < dj[j] && bi[i] != 0.0)
			{
				lowVal = bi[i];
				dj[j] -= bi[i];
				bi[i] = 0.0;

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

				// Test on Degeneration
				if (bi[i] == 0.0 && dj[j] == 0.0 && i < m - 1 && j < n - 1)
				{
					isBasic[i][j+1] = 1;
					basicsEnd->val = lowVal;
					basicsEnd->i = i;
					basicsEnd->j = j+1;

					basicsEnd->nextCus = cusBasics[j+1];
					if (cusBasics[j+1] != NULL) cusBasics[j+1]->prevCus = basicsEnd;
					basicsEnd->nextFac = facBasics[i];
					if (facBasics[i] != NULL) facBasics[i]->prevFac = basicsEnd;

					facBasics[i] = basicsEnd;
					basicsEnd->prevCus = NULL;
					cusBasics[j+1] = basicsEnd;
					basicsEnd->prevFac = NULL;

					basicsEnd++;

					i++;
					j++;
				}
				else
				{
					// Skip to next Customer
					i++;
				}

			}
		}
	}

	/**********************
	LCM initialization method
	**********************/
	void _initLCM(
		double *bi,
		double *dj,
		TsBasic *basicsEnd,
		TsBasic **facBasics,
		TsBasic **cusBasics,
		bool ** isBasic,
		int m,
		int n,
		ublas::matrix<unsigned int> &cij_column_sorted,
		int *facilities)
	{
		unsigned int i = 0, j = 0;
		unsigned int indx_i = 0;
		double lowVal = 0.0;

		while (j < cij_column_sorted.size2())
		{
			indx_i = cij_column_sorted(i, j);

			// Can shift Capacity to fullfill Customer Demand
			if (bi[indx_i] != 0.0 && dj[j] != 0.0)
			{
				// More capacity than Demand
				if (bi[i] >= dj[j] && dj[j] != 0.0)
				{
					lowVal = dj[j];
					bi[i] -= dj[j];
					dj[j] = 0.0;

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

					// Skip to next Customer
					j++;
					i = 0;
				}
				// Less Capacity than Demand
				else if (bi[i] < dj[j] && bi[i] != 0.0)
				{
					lowVal = bi[i];
					dj[j] -= bi[i];
					bi[i] = 0.0;

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

					// Skip to next Customer
					i++;
				}
			}
		}

		// If its in the fictional Column
		if (n > cij_column_sorted.size2() && j == cij_column_sorted.size2())
		{
			for (i = 0; i != cij_column_sorted.size1(); ++i)
			{
				// If Capacity is not empty ship to customer
				if (bi[i] != 0.0)
				{
					lowVal = bi[i];
					bi[i] = 0.0;
					
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
			}
		}
	}
}