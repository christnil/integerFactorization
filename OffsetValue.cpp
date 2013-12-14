/*
 * OffsetValue.cpp
 *
 *  Created on: Oct 25, 2012
 *      Author: christoffer
 */

#include "OffsetValue.h"

OffsetValue::OffsetValue(int poffset, int pvalue, unsigned char plogvalue) {
	offset = poffset;
	value = pvalue;
	while (offset < 0)
		offset += value;
	logvalue = plogvalue;

}

OffsetValue::~OffsetValue() {

}

