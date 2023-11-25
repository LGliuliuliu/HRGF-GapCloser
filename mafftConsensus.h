#ifndef MAFFTCONSENSUS_H_INCLUDED 
#define MAFFTCONSENSUS_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "scaffoldSet.h"


typedef struct AlignmentArray {
	int* weight;
}AlignmentArray;

typedef struct AlignmentArrayHead {
	AlignmentArray* alignmentArray;

}AlignmentArrayHead;




bool GetMafftConsensusSequenceForOneGap(ScaffoldSetHead* scaffoldSetHead, long int i, long int j, char* resultOutPutDirectory);

bool GetConsensusSequence(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory);

void MAFFT(char* tempSequence, char* alignmentFile);

ContigSetHead* GetContigSetForAlignmentFile(char* alignmentFile);

long int VoteMajority(int* weight, int count);

char* GetConsensusSequenceFromOne(ContigSetHead* contigSetHead);

void ReverseSequenceReadToGap(ReadToGap* readToGap, long int count);

bool GetMafftLeftAndRightSequenceTarget(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory);



#endif
