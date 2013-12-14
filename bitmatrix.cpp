#include "bitmatrix.h"

bitmatrix::bitmatrix(void)
{
}

bitmatrix::bitmatrix(int p_rows, int p_cols)
{
	rows = p_rows;
	cols = p_cols;
	colwords = p_cols / WORDLENGTH;
	if (colwords == (p_cols - 1) / WORDLENGTH)
		colwords++;
	m = new unsigned int[rows * colwords];
	leadingZero = new int[cols];
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < colwords; ++j) {
			*(m + (i * colwords) + j) = 0;
		}
	}
}

void bitmatrix::setbit_1(int row, int col)
{
	int wordinrow = col / WORDLENGTH;
	unsigned int mask = 1 << (WORDLENGTH - 1 - (col % WORDLENGTH));
	*(m + (row * colwords) + wordinrow) |= mask;
}
void bitmatrix::setbit_0(int row, int col)
{
	int wordinrow = col / WORDLENGTH;
	unsigned int mask = 1 << (WORDLENGTH - 1 - (col % WORDLENGTH));
	*(m + (row * colwords) + wordinrow) &= ~mask;
}
void bitmatrix::setcolumn(int colindex, std::vector<int> values)
{
	if (values.size() != rows) throw(1);
	int wordinrow = colindex / WORDLENGTH;
	unsigned int setmask = 1 << (WORDLENGTH - 1 - (colindex % WORDLENGTH));
	unsigned int unsetmask = ~setmask;
	for (int row = 0; row < rows; ++row)
	{
		if (values[row] % 2 == 1)
			*(m + (row * colwords) + wordinrow) |= setmask;
		else
			*(m + (row * colwords) + wordinrow) &= unsetmask;
	}
}
void bitmatrix::setcolumn(int colindex, int * values)
{
	int wordinrow = colindex / WORDLENGTH;
	unsigned int setmask = 1 << (WORDLENGTH - 1 - (colindex % WORDLENGTH));
	unsigned int unsetmask = ~setmask;
	for (int row = 0; row < rows; ++row)
	{
		if (values[row] % 2 == 1)
			m[ + (row * colwords) + wordinrow] |= setmask;
		else
			m[ + (row * colwords) + wordinrow] &= unsetmask;
	}
}

char bitmatrix::getbitchar(int row, int col) const
{
	int wordinrow = col / WORDLENGTH;
	unsigned int mask = 1 << (WORDLENGTH - 1 - (col % WORDLENGTH));
	if (*(m + (row * colwords) + wordinrow) & mask)
		return '1';
	return '0';
}
bool bitmatrix::getBit(int row, int col)
{
	int wordinrow = col / WORDLENGTH;
	unsigned int mask = 1 << (WORDLENGTH - 1 - (col % WORDLENGTH));
	if (*(m + (row * colwords) + wordinrow) & mask)
		return true;
	return false;
}

void bitmatrix::switchrows(int row1, int row2)
{
	if (row1 == row2) return;
	unsigned int tmp;
	for (int i = 0; i < colwords; ++i) {
		tmp = *(m + (row1 * colwords) + i);
		*(m + (row1 * colwords) + i) = *(m + (row2 * colwords) + i);
		*(m + (row2 * colwords) + i) = tmp;
	}
	/*std::cout << "switching rows: " << row1 << " : " << row2 << std::endl;
	print();
	std::cout << std::endl;*/
}

void bitmatrix::xor_rows(int row, int targetrow)
{
	for (int i = 0; i < colwords; ++i) {
		*(m + (targetrow * colwords) + i) ^= *(m + (row * colwords) + i);
	}
	/*std::cout << "xor rows: " << row << " : " << targetrow << std::endl;
	print();
	std::cout << std::endl;*/
}

void bitmatrix::gauss()
{
	int firstfreerow = 0;
	int row;
	for (int col = 0; col < cols; ++col)
	{
		for (row = firstfreerow; row < rows; ++row)
		{
			if (getBit(row, col))
			{
				switchrows(firstfreerow, row);
				break;
			}
		}
		if (getBit(firstfreerow, col))
		{
			for (++row; row < rows; ++row)
			{
				if (getBit(row, col))
				{
					xor_rows(firstfreerow, row);
				}
			}
			++firstfreerow;
		}
	}
	//set free vector
	for (int i = 0; i < cols; i++)
	{
		leadingZero[i] = -1;
	}
	numberOfSolutions = cols;
	int row1 = 0;
	int col = 0;
	for ( ; row1 < rows && col < cols; )
	{
		while (col < cols && !getBit(row1, col))
		{
			++col;
		}
		if (col < cols && getBit(row1, col))
		{
			leadingZero[col] = row1;
			--numberOfSolutions;
		}
		++row1;
		++col;
	}
	if (numberOfSolutions > 30)
		numberOfSolutions = 30;
	numberOfSolutions = (1 << numberOfSolutions);
	currentsolution = 0;
}

std::vector<int> bitmatrix::getNextSolution(void)
{
	currentsolution++;
	int bitmask = 1;
	unsigned int * tmp = new unsigned int[colwords];
	for (int i = 0; i < colwords; i++)
	{
		tmp[i] = 0;
	}
	std::vector<int> a(cols, 0);
	for (int col = cols - 1; col >= 0; --col)
	{
		if (leadingZero[col] == -1)
		{
			if (currentsolution & bitmask)
			{
				a[col] = 1;
				tmp[col / WORDLENGTH] |= 1 << (WORDLENGTH - 1 - (col % WORDLENGTH));
			}
			else
			{
				a[col] = 0;
			}
			bitmask = bitmask << 1;
		}
		else
		{
			if (oddBits(m + (colwords * leadingZero[col]), tmp, colwords))
			{
				a[col] = 1;
				tmp[col / WORDLENGTH] |= 1 << (WORDLENGTH - 1 - (col % WORDLENGTH));
			}
		}
	}
	delete tmp;
	return a;
}

bool bitmatrix::oddBits(unsigned int * firstWord, unsigned int * mask, int numWords)
{
	bool odd = false;
	unsigned int tmp;
	for (int i = 0; i < numWords; i++)
	{
		tmp = firstWord[i] & mask[i];
		for (long int j = 1; j != 0; j = j << 1)
		{
			if ((tmp & j) != 0)
				odd = !odd;
		}
	}
	return odd;
}

bool bitmatrix::hasMoreSolutions(void)
{
	return currentsolution < (numberOfSolutions - 1);
}

std::ostream &operator<<(std::ostream &out, const bitmatrix &b)
{
	for (int i = 0; i < b.rows; ++i) {
		for (int j = 0; j < b.cols; ++j) {
			out << b.getbitchar(i, j) << " ";
		}
		out << std::endl;
	}
	out << "free vector:" << std::endl;
	for (int i = 0; i < b.cols; i++)
	{
		out << b.leadingZero[i] << " ";
		/*if (b.leadingZero[i] == -1)
			out << "1 ";
		else 
			out << "0 ";*/
	}
	out << std::endl;
	return out;
}

void bitmatrix::print()
{
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			std::cout << getbitchar(i, j) << " ";
		}
		std::cout << std::endl;
	}
}

bitmatrix::~bitmatrix(void)
{
	delete[] m;
	delete[] leadingZero;
}
