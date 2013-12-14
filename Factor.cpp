//============================================================================
// Name        : Factor.cpp
// Author      : Christoffer Nilsson
// Version     :
// Copyright   : 
// Description : Prime factoring C++
//============================================================================

#include <iostream>
#include <vector>
#include <stack>
#include "Sieve.h"

#define NUMBERS 100

using namespace std;

bool factor(mpz_class, vector<mpz_class>&);
bool perfectpower(mpz_class, vector<mpz_class>&);

bool perfectpower(mpz_class n, vector<mpz_class>& factors)
{
	if (n == 1) return true;
	int power;
	mpz_class tmp;
	for (power = 8; power > 1; --power)
	{
		mpz_root(tmp.get_mpz_t(), n.get_mpz_t(), power);
		mpz_pow_ui(tmp.get_mpz_t(), tmp.get_mpz_t(), power);
		if (tmp == n) break;
	}
	mpz_root(tmp.get_mpz_t(), n.get_mpz_t(), power);
	if (mpz_probab_prime_p(tmp.get_mpz_t(), 10))
	{
		for (int i = 0; i < power; i++)
			factors.push_back(tmp);
		return true;
	}
	else
	{
		vector<mpz_class> tmpvec;
		if (!factor(tmp, tmpvec))
			return false;
		for (vector<mpz_class>::iterator it = tmpvec.begin();
				it != tmpvec.end(); ++it)
		{
			for (int i = 0; i < power; ++i)
				factors.push_back(*it);
		}
		return true;
	}
}

bool factor(mpz_class n, vector<mpz_class>& factors)
{
	if (n == 1) return true;
	if (mpz_probab_prime_p(n.get_mpz_t(), 10))
	{
		factors.push_back(n);
		return true;
	}
	if (mpz_perfect_power_p(n.get_mpz_t()))
	{
		return perfectpower(n, factors);
	}
	Sieve s(n);
	vector<mpz_class> factorsfound = s.Factor();
	vector<mpz_class> nonprimes;
	if (factorsfound.size() == 0)
		return false;
	bool foundOne = false;
	for (vector<mpz_class>::iterator it = factorsfound.begin();
			it != factorsfound.end(); ++it)
	{
		if (mpz_probab_prime_p(it->get_mpz_t(), 10))
		{
			factors.push_back(*it);
			n = n / *it;
			foundOne = true;
		}
	}
	if (foundOne)
		return factor(n, factors);
	else
	{
		if (!factor(*factorsfound.begin(), factors))
			return false;
		if (!factor(n / *factorsfound.begin(), factors))
			return false;
	}
	return true;
}

int main() {
	mpz_class tofactor[NUMBERS];
	vector<mpz_class> factors[NUMBERS];
	bool done[NUMBERS];
	string s;
	for (int i = 0; i < NUMBERS; ++i)
	{
		cin >> s;
		mpz_class input(s);
		for (int j = 0; j < lowprimes::n; ++j)
		{
			if (input % lowprimes::primearray[j] == 0)
			{
				input = input / lowprimes::primearray[j];
				factors[i].push_back(lowprimes::primearray[j]);
				--j;
			}
		}
		if (input == 1)
		{
			done[i] = true;
		}
		else
		{
			done[i] = false;
			tofactor[i] = input;
		}
	}
	//factor with QS
	mpz_class tmp;
	for (int i = 0; i < NUMBERS; ++i)
	{
		if (done[i]) continue;
		done[i] = factor(tofactor[i], factors[i]);
	}
	//print result
	for (int i = 0; i < NUMBERS; ++i)
	{
		if (done[i])
		{
			for (vector<mpz_class>::iterator it = factors[i].begin();
					it != factors[i].end(); ++it)
			{
				cout << *it << endl;
			}
		}
		else
		{
			cout << "fail" << endl;
		}
		cout << endl;
	}
	return 0;
}
