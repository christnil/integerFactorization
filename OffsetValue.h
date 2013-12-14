/*
 * OffsetValue.h
 *
 *  Created on: Oct 25, 2012
 *      Author: christoffer
 */

#ifndef OFFSETVALUE_H_
#define OFFSETVALUE_H_
#include <gmp.h>
#include <gmpxx.h>

struct OffsetValue {
	OffsetValue(int poffset, int pvalue, unsigned char plogvalue);
	virtual ~OffsetValue();
	int offset;
	int value;
	unsigned char logvalue;
};

#endif /* OFFSETVALUE_H_ */
