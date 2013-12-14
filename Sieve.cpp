/*
 * Sieve.cpp
 *
 *  Created on: Oct 25, 2012
 *      Author: christoffer
 */

#include "Sieve.h"

Sieve::Sieve(mpz_class pN)
{
	N = pN;
	mpz_root(N_sqrt.get_mpz_t(), N.get_mpz_t(), 2);
	floatingthreshhold = THRESHHOLD;
	numFoundInSweepp = 0;
	numtriedInSweep = 0;
	SetUp();
}

Sieve::~Sieve() {
	delete matrixptr;
	delete[] array;
}

void Sieve::updateThreshhold()
{
	float sRate = ((float)numFoundInSweepp) / ((float)numtriedInSweep);
	numFoundInSweepp = 0;
	numtriedInSweep = 0;
	if (sRate > 0.9)
	{
		floatingthreshhold += 0.5;
		return;
	}
	if (sRate < 0.8)
	{
		floatingthreshhold -= ((0.8 - sRate) * 1.5);
	}
}

void Sieve::SetUp(void)
{
#ifdef PRINT
	std::cout << "calculating factor base bound" << std::endl;
	CalculateB();
	std::cout << "bound = " << B << std::endl;
	std::cout << "creating factor base" << std::endl;
	CreateFactorBase();
	std::cout << "factor base length = " << factorbase.size() << std::endl;
	std::cout << "calculating offsets" << std::endl;
	CalculateOffsets2();
	std::cout << offsets.size() << " offsets found" << std::endl;
#endif
#ifndef PRINT
	CalculateB();
	CreateFactorBase();
	CalculateOffsets2();
#endif
	numsmooth = 0;
	numtried = 0;
	numsmoothrequired = factorbase.size() + NUMEXTRA;
	matrixptr = new bitmatrix(factorbase.size(), numsmoothrequired);

}

std::vector<mpz_class> Sieve::Factor()
{
	unsigned long int first = 0;
	int iter = 0;

#ifdef PRINT
					std::cout << "beginning to factor"<< std::endl;
					std::cout << "num smooth required: " << numsmoothrequired << std::endl;
					std::cout << "MAXSWEEPS: " << MAXSWEEPS << std::endl;
#endif
	while (numsmooth < numsmoothrequired && iter++ < MAXSWEEPS)
	{
		Sweep(first);
		updateThreshhold();
		first += SIEVESIZE;
	}
#ifdef PRINT
					std::cout << iter <<  " sweeps done"<< std::endl;
#endif
	std::vector<mpz_class> factors;
	std::vector<int> nullvec;
#ifdef PRINT
					std::cout << "num smooth found: " << numsmooth << std::endl;
					std::cout << "num tried: " << numtried << std::endl;
#endif
	if (numsmooth == numsmoothrequired)
	{
		mpz_class a, b, ab, factor1;
		matrixptr->gauss();
		for (int i = 0; i < 20 && matrixptr->hasMoreSolutions(); ++i)
		{
			nullvec = matrixptr->getNextSolution();
			a = 1;
			b = 1;
			for (int j = 0; j < nullvec.size(); j++)
			{
				if (nullvec[j] == 1)
				{
					a = (a * smooth_x[j]);
					b = (b * smooth_y[j]);
				}
			}
			mpz_root(b.get_mpz_t(), b.get_mpz_t(), 2);
			ab = a - b;
			mpz_abs(ab.get_mpz_t(), ab.get_mpz_t());
			mpz_gcd(factor1.get_mpz_t(), ab.get_mpz_t(), N.get_mpz_t());

			if (factor1 != 1 && factor1 != N)
			{
				bool exists = false;
				for (std::vector<mpz_class>::iterator it = factors.begin(); it != factors.end(); ++it)
				{
					if ((*it) == factor1)
					{
						exists = true;
						break;
					}
				}
				if (!exists)
				{
					factors.push_back(factor1);
					i = 0;
				}
			}
		}
	}
	return factors;
}

void Sieve::Sweep(unsigned long int start_X)
{
	unsigned long int last = start_X + SIEVESIZE - 1;
	mpz_class qlast = Q(last);
	unsigned char logval = std::log(qlast.get_d()) + 0.9999;
	for (int i = 0; i < SIEVESIZE; ++i)
		sieve[i] = logval;

	for (std::vector<OffsetValue>::iterator it = offsets.begin();
			it != offsets.end(); ++it)
	{
		while (it->offset < SIEVESIZE)
		{
			sieve[it->offset] -= it->logvalue;
			it->offset += it->value;
		}
		it->offset -= SIEVESIZE;
	}

	for (int i = 0; i < SIEVESIZE && numsmooth < numsmoothrequired; ++i)
	{
		if (sieve[i] < floatingthreshhold)
			TryFactor(start_X + i);
	}
}

void Sieve::TryFactor(unsigned long int x)
{
	if (numsmooth == numsmoothrequired) return;
	numtried++;
	numtriedInSweep++;
	mpz_class mpzx = N_sqrt + x;
	mpz_class mpzy = Q(x);
	mpz_class copy_y(mpzy);

#ifdef PRINT
					//std::cout << "trying to factor: " << mpzy << std::endl;
#endif
	for (int i = 0; i < factorbase.size(); ++i)
		array[i] = 0;

	int i = 0;
	for (std::vector<unsigned long int>::iterator it = factorbase.begin();
			it != factorbase.end(); ++it, ++i)
	{
		while (copy_y % *it == 0)
		{
			copy_y = copy_y / *it;
			array[i]++;
		}
	}

	if (copy_y == 1)
	{
#ifdef PRINTALL
					std::cout << "found smooth: " << mpzy << " xval: " << mpzx <<
							" index: " << x << std::endl;
#endif
		matrixptr->setcolumn(numsmooth, array);
		smooth_x.push_back(mpzx);
		smooth_y.push_back(mpzy);
		numsmooth++;
		numFoundInSweepp++;
	}
}

void Sieve::CalculateB()
{
	double n_double = N.get_d();
	B = C * std::exp(0.5 * std::sqrt(std::log(n_double) * std::log(std::log(n_double)))) + C2;
}

void Sieve::CreateFactorBase(void)
{
	factorbase.push_back(2);
	mpz_class tmpmpz;
	for (int i = 1; i < lowprimes::n && lowprimes::primearray[i] < B; i++)
	{
		tmpmpz = lowprimes::primearray[i];
		if (mpz_legendre(N.get_mpz_t(), tmpmpz.get_mpz_t()) == 1)
		{
#ifdef PRINTALL
			std::cout << lowprimes::primearray[i] << " ";
#endif
			factorbase.push_back(lowprimes::primearray[i]);
		}
	}
#ifdef PRINTALL
	std::cout << std::endl;
#endif

	array = new int[factorbase.size()];
}

void Sieve::CalculateOffsets(void)
{
	unsigned long int tmpSize = SIEVESIZE;// std::min(1000, SIEVESIZE);
	mpz_class * tmpArray = new mpz_class[tmpSize];
	for (unsigned long int i = 0; i < tmpSize; ++i)
	{
		tmpArray[i] = Q(i);
	}

	unsigned long int ack;
	float logval;
	for (std::vector<unsigned long int>::iterator it = factorbase.begin();
			it != factorbase.end(); ++it)
	{
		logval = std::log((float)*it);
		for (ack = *it; ack < B; ack *= *it)
		{
#ifdef PRINTALL
					std::cout << "trying ack: " << ack << std::endl;
#endif
			for (unsigned long int i = 0; i < tmpSize; ++i)
			{
				if (tmpArray[i] % ack == 0)
				{
#ifdef PRINTALL
					std::cout << "value: " << ack << " offset: " << i << " approx: " <<
							logval << std::endl;
#endif
					OffsetValue offv(i, ack, logval);
					offsets.push_back(offv);
					for (unsigned long int j = i + 1; j < tmpSize && j < i + ack; j++)
					{
						if (tmpArray[j] % ack == 0)
						{
#ifdef PRINTALL
					std::cout << "value: " << ack << " offset: " << j << " approx: " <<
							logval << std::endl;
#endif
							OffsetValue offv2(j, ack, logval);
							offsets.push_back(offv2);
							if (it != factorbase.begin())
								break;
						}
					}
					break;
				}
			}
		}
	}

	delete[] tmpArray;
}



#define mpz_rshift(A,B,l) mpz_tdiv_q_2exp(A, B, l)

int root(mpz_class& result, const mpz_class arg, const mpz_class prime)
{
	mpz_class y, b, t;
	unsigned int r, m;

	if (arg % prime == 0) {
		result = 0;
		return 1;
	}
	if (mpz_legendre(arg.get_mpz_t(), prime.get_mpz_t()) == -1)
		return -1;

	b = 0;
	t = 0;
	y = 2;
	while(mpz_legendre(y.get_mpz_t(), prime.get_mpz_t()) != -1)
	{
		y++;
	}

	result = prime - 1;
	r = mpz_scan1(result.get_mpz_t(), 0);
	mpz_rshift(result.get_mpz_t(), result.get_mpz_t(), r);

	mpz_powm(y.get_mpz_t(), y.get_mpz_t(), result.get_mpz_t(), prime.get_mpz_t());

	mpz_rshift(result.get_mpz_t(), result.get_mpz_t(), 1);
	mpz_powm(b.get_mpz_t(), arg.get_mpz_t(), result.get_mpz_t(), prime.get_mpz_t());

	result = (arg * b) % prime;

	b = (result * b) % prime;

	while(b != 1) {
		t = (b * b) % prime;
		for(m = 1; t != 1; m++) {
			t = (t * t) % prime;
		}

		t = 0;
		mpz_setbit(t.get_mpz_t(), r - m - 1);
		mpz_powm(t.get_mpz_t(), y.get_mpz_t(), t.get_mpz_t(), prime.get_mpz_t());

		y = t * t;

		r = m;

		result = (result * t) % prime;

		b = (b * y) % prime;
	}
	return 1;
}

mpz_class Sieve::Tonelli_Shanks(mpz_class p, mpz_class n)
{
	mpz_class R;
	int res = root(R, n, p);
	//std::cout << "retval: " << res << std::endl;
	return R;
}
void Sieve::CalculateOffsets2(void)
{
	if (factorbase.size() == 0) return;
	if (Q(0) % 2 == 0)
	{
		OffsetValue off2(0, 2, std::log(2.0));
		offsets.push_back(off2);
		CalculateOffsetsRec(0, 2, 2, std::log(2.0));
	}
	else
	{
		OffsetValue off2(1, 2, std::log(2.0));
		offsets.push_back(off2);
		CalculateOffsetsRec(1, 2, 2, std::log(2.0));
	}
	for (std::vector<unsigned long int>::iterator it = factorbase.begin() + 1;
			it != factorbase.end(); ++it)
	{
		mpz_class retv = Tonelli_Shanks(*it, N);
		mpz_class cof = *it - retv;
		retv = (retv - N_sqrt) % *it;
		if (retv < 0)
			retv += *it;
		cof = (cof - N_sqrt) % *it;
		if (cof < 0)
			cof += *it;
		if (retv < SIEVESIZE && retv >= 0)
		{
			OffsetValue off(retv.get_ui(), *it, std::log((float)*it));
			offsets.push_back(off);
			CalculateOffsetsRec(retv.get_ui(), *it, *it, std::log((float)*it));
		}
		if (cof < SIEVESIZE && cof >= 0)
		{
			OffsetValue off(cof.get_ui(), *it, std::log((float)*it));
			offsets.push_back(off);
			CalculateOffsetsRec(cof.get_ui(), *it, *it, std::log((float)*it));
		}
	}
}
void Sieve::CalculateOffsetsRec(int prevExpOff, int step, int factor, float logvalue)
{
	step = step * factor;
	if (step > B) return;
	if (prevExpOff > -1)
	{
		for (int i = prevExpOff, j = 0; i < SIEVESIZE && j < factor; i += step, j++)
		{
			if (Q(i) % step == 0)
			{
				OffsetValue off(i, step, logvalue);
				offsets.push_back(off);
				CalculateOffsetsRec(i, step, factor, logvalue);
				if (factor != 2)
					break;
			}
		}
	}
}

mpz_class Sieve::Q(unsigned long int  x)
{
	mpz_class tmp = N_sqrt + x;
	tmp = tmp * tmp;
	tmp = tmp - N;
	return tmp;
}

