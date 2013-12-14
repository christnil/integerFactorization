#pragma once
#include <iostream>
#include <vector>

#define WORDLENGTH 32

class bitmatrix
{
public:
	friend std::ostream& operator<<(std::ostream& os, const bitmatrix& dt);
	
	bitmatrix(void);
	bitmatrix(int , int);
	~bitmatrix(void);
	
	void setbit_1(int, int);
	void setcolumn(int colindex, std::vector<int> values);
	void setcolumn(int colindex, int * values);
	void setbit_0(int, int);
	void switchrows(int row1, int row2);
	void xor_rows(int row, int targetrow);
	void gauss(void);
	void print(void);
	
	char getbitchar(int, int) const;
	bool getBit(int row, int col);

	std::vector<int> getNextSolution(void);
	bool hasMoreSolutions(void);
	bool oddBits(unsigned int * firstWord, unsigned int * mask, int numWords);

private:
	int * leadingZero;
	int currentsolution;
	int numberOfSolutions;
	int rows;
	int cols;
	unsigned int * m;
	int colwords;
};

