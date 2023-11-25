#ifndef SCAFFOLDSET_H_INCLUDED 
#define SCAFFOLDSET_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "readSet.h"
#include"contigSet.h"

typedef struct ReadToGap {

	char* read;
	long int readLength;
	long int weigth;
	long int leftContigDiff;
	long int rightContigDiff;
	long int leftContigOverlap;
	long int rightContigOverlap;

}ReadToGap;


typedef struct HifiGapSet {
	ReadToGap* leftReadSet;
	ReadToGap* rightReadSet;
	ReadToGap* spanReadSet;
	long int leftReadCount;
	long int rightReadCount;
	long int spanReadCount;
	long int leftTooLengthCount;
	long int rightTooLengthCount;
	long int spanFlag;
	char* spanConsensusSequence;
	char* leftConsensusSequence;
	char* rightConsensusSequence;
}HifiGapSet;

typedef struct GapSet {
	long int gapStartCoordinate;
	long int gapEndCoordinate;
	long int gapLength;
	long int estimatedGapDistance;
	long int estimatedReadCount;
}GapSet;

typedef struct ContigPosition {
	long int contigStartCoordinate;
	long int contigEndCoordinate;
	long int contigLength;
}ContigPosition;
typedef struct ScaffoldSet {
	char* scaffoldName;
	char* scaffold;
	long int scaffoldLength;
	GapSet* gapSet;
	ContigPosition* contigPostion;
	HifiGapSet* hifiGapSet;
	long int gapCount;
	long int contigCount;

}ScaffoldSet;

typedef struct ScaffoldSetHead {
	ScaffoldSet* scaffoldSet;
	long int scaffoldCount;
	long int allGapCount;
}ScaffoldSetHead;

typedef struct ReadAlign {
	long int readIndex;
	long int readLength;
	long int readStartCoordinate;
	long int readEndCoordinate;
	long int orinetation;
	long int minOverlap;
	double identityRatio;
}ReadAlign;

typedef struct GapReadAlignSet {
	long int spanReadAlignCount;
	ReadAlign* spanReadAlign;
	long int leftReadAlignCount;
	ReadAlign* leftReadAlign;
	long int rightReadAlignCount;
	ReadAlign* rightReadAlign;
}GapReadAlignSet;

typedef struct Scaffold {
	GapReadAlignSet* gapReadAlignSet;
	long int gapCount;
}Scaffold;

typedef struct ScaffoldHead {
	Scaffold* scaffold;
	long int scaffoldCount;

}ScaffoldHead;



void DelectShortContig(char* scaffoldSetFile, char* longScaffoldSetFile);

ScaffoldSetHead* GetScaffoldSetFromScaffoldFile(char* scaffoldSetFile);

GapSet* GetGapSetFromOneScaffold(char* scaffold, long int& gapCount);

void GetContigSetFromScaffoldSetHead(ScaffoldSetHead* scaffoldSetHead);

void OutPutContigSetInScaffoldSetHead(ScaffoldSetHead *scaffoldSetHead,char *contigSetFile);

void GetEstimatedGapLength(ScaffoldSetHead* scaffoldSetHead, ContigSetHead* contigSetHead, char* file, char* line, int maxSize);

void GetHiFiInGapFromReadPoinst(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, char* file, char* line, int maxSize, FILE* hifiFileFP);

long int GetGapIndexFromContigIndex(ScaffoldSetHead* scaffoldSetHead, long int leftContigIndex, long int rightContigIndex, long int& scaffoldIndex, long int& gapIndex);

void GetReadPoinstionInGap(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, char* file, char* line, int maxSize, FILE* longReadFileFP);

long int GetGapIndexFromContigIndexOne(ScaffoldSetHead* scaffoldSetHead, long int contigIndex, bool right, long int& scaffoldIndex, long int& gapIndex);

void GetHiFiInGap(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, char* longReadFile, char* positionFile, char* line, int maxSize, char* outputAddress);

long int GetContigIndexFromGapIndex(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, long int& gapLeftContigindex, long int& gapRightContigIndex);

void GetSupplementReadPoinstionInGap(ScaffoldSetHead* scaffoldSetHead, ContigSetHead* contigSetHead, char* file, char* line, int maxSize, FILE* hifiFileFP);

void GetSupplementHiFiInGap(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, char* hifiFile, char* positionFile, char* line, int maxSize, char* outputAddress);







#endif