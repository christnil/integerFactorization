/*
 * Sieve.h
 *
 *  Created on: Oct 25, 2012
 *      Author: christoffer
 */

#ifndef SIEVE_H_
#define SIEVE_H_

#include "lowprime.h"

#include "bitmatrix.h"
#include "OffsetValue.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>

#define SIEVESIZE 25000
#define THRESHHOLD 12
#define MAXSWEEPS 60000
#define C 3
#define C2 500
#define NUMEXTRA 30

//#define PRINT


class Sieve {
public:
	Sieve(mpz_class pN);
	virtual ~Sieve();
	std::vector<mpz_class> Factor();

private:
	//functions
	mpz_class Tonelli_Shanks(mpz_class p, mpz_class n);
	void SetUp(void);
	void CalculateB(void);
	void CreateFactorBase(void);
	void CalculateOffsets(void);
	void CalculateOffsets2(void);
	void CalculateOffsetsRec(int prevExpOff, int step, int factor, float logvalue);


	void Sweep(unsigned long int  start_X);
	void TryFactor(unsigned long int x);
	mpz_class Q(unsigned long int  x);


	// class members
	mpz_class N, N_sqrt;
	unsigned long int B;
	unsigned char sieve[SIEVESIZE];
	std::vector<unsigned long int> factorbase;
	std::vector<OffsetValue> offsets;
	bitmatrix * matrixptr;

	int numsmooth;
	int numsmoothrequired;
	int numtried;
	std::vector<mpz_class> smooth_x;
	std::vector<mpz_class> smooth_y;
	int * array;

	void updateThreshhold(void);
	float floatingthreshhold;
	int numtriedInSweep, numFoundInSweepp;
};

#endif /* SIEVE_H_ */
