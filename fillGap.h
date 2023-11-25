#ifndef FILLGAP_H_INCLUDED 
#define FILLGAP_H_INCLUDED 


#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "readSet.h"
#include "contigSet.h"
#include "scaffoldSet.h"
#include "mafftConsensus.h"



void FillingGap(ScaffoldSetHead* scaffoldSetHead, char* resultOutPutDirectory);

char* CombineLeftAndRightConsensusSequence(char* leftConsensusSequence, char* rightConsensusSequence, long int gapDistance);

void OutputFillingResultFromScaffoldSetHead(ScaffoldSetHead* scaffoldSetHead, char* result);

void OutputFillingResultFromScaffoldSetHead1(ScaffoldSetHead* scaffoldSetHead, char* fillingResultFile);

void TransferToLower(char* contig);

long int GetContigIndexFromGapIndex(ScaffoldSetHead* scaffoldSetHead, long int i, long int j);

#endif