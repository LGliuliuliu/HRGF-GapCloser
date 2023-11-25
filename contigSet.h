#ifndef CONTIG_H_INCLUDED 
#define CONTIG_H_INCLUDED 


#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

using namespace std;

typedef struct ContigSet {
	char* contig;
	char* contigName;
	long int contigLength;
}ContigSet;

typedef struct ContigSetHead {
	ContigSet* contigSet;
	long int contigCount;
}ContigSetHead;

ContigSetHead* GetContigSetFromContigSetFile(char* contigSetFile);

char* ReverseComplement(char* temp);

long int GetMinValue(long int a, long int b);

void ReverseSequence(char* temp);




#endif