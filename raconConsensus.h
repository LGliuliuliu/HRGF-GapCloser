#ifndef RACONCONSENSUS_H_INCLUDED 
#define RACONCONSENSUS_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "contigSet.h"
#include "readSet.h"
#include "scaffoldSet.h"
#include"mafftConsensus.h"


typedef struct ReadSortCountForSimilarSet {
	long int readIndex;
	long int flag;
}ReadSortCountForSimilarSet;

typedef struct ReadSortForSimilar {
	long int readSortForSimilarCount;
	ReadSortCountForSimilarSet* readSortCountForSimilarSet;
}ReadSortForSimilar;

typedef struct ReadSortCountForLengthSet {
	long int readIndex;
	long int flag;
}ReadSortCountForLengthSet;

typedef struct ReadSortForLength {
	ReadSortCountForLengthSet* readSortCountForLengthSet;
	long int readSortForLengthCount;
	long int aveReadLength;
	ReadSortForSimilar* readSortCountForSimilar;
	long int sortSimilarCount;

}ReadSortForLength;

typedef struct ReadSort {
	ReadSortForLength* readSortForLength;
	long int sortLengthCount;
}ReadSort;


using namespace std;

bool GetRaconConsensusSequenceForOneGap(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, char* resultOutPutDirectory);

bool GetSpanSequenceTarget(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex,long int gapIndex, int type,char* resultOutPutDirectory);

bool GetLeftAndRightSequenceTarget(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, char* dataBase, char* overlapFile, char* targetFile, char* sequenceFile, char* consensusFile, int type);

char* GetConsensusSequenceFromFile(char* file);

bool RunRacon(char* overlapFile, char* targetFile, char* sequenceFile, char* consensusFile);

bool GetleftAndRightConsenSequence(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory);

bool GetConsenSequence(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory);

#endif 
