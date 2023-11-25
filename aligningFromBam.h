#ifndef aligningFromBam_H_INCLUDED 
#define aligningFromBam_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "contigSet.h"
#include "scaffoldSet.h"

using namespace std;
using namespace BamTools;

typedef struct AligningResult {
	int readStartPosition;
	int readEndPosition;
	int contigStartPosition;
	int contigEndPosition;
	int contigIndex;
	long int readIndex;
	int overlapLength;
	bool orientation;
}AligningResult;

typedef struct AligningResultHead {
	AligningResult* aligningResult;
	int allocateAligningResultCount;
	int aligningResultCount;
	int aligningShortContigResultCount;
}AligningResultHead;



int GetAligningResult(ScaffoldSetHead* scaffoldSetHead,ContigSetHead* contigSetHead, char* aligningResultBam, char* file);

bool GetAligningResultOneLine(AligningResultHead* aligningResultHead, BamAlignment alignment, ContigSetHead* contigSetHead, long int index);

void OutPutAligningResultOneLine(ScaffoldSetHead* scaffoldSetHead,AligningResultHead* aligningResultHead, ContigSetHead* contigSetHead, FILE* fp, long int readIndex, long int readLength);

int GetAligningResult1(ContigSetHead* contigSetHead, char* aligningResultBam, char* file);

bool GetAligningResultOneLine1(AligningResultHead* aligningResultHead, BamAlignment alignment, ContigSetHead* contigSetHead, long int index);

void OutPutAligningResultOneLine1(AligningResultHead* aligningResultHead, ContigSetHead* contigSetHead, FILE* fp, long int readIndex, long int readLength);


#endif
