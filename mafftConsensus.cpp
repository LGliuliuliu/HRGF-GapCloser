#ifndef MAFFTCONSENSUS_CPP_INCLUDED 
#define MAFFTCONSENSUS_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include"mafftConsensus.h"
#include "scaffoldSet.h"

using namespace std;

bool GetMafftConsensusSequenceForOneGap(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex,char* resultOutPutDirectory) {

	if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanReadCount > 0) {
		GetConsensusSequence(scaffoldSetHead, scaffoldIndex, gapIndex, 0, resultOutPutDirectory);
		return true;
	}
	else {
		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftReadCount > 0) {
			GetMafftLeftAndRightSequenceTarget(scaffoldSetHead, scaffoldIndex, gapIndex, 1, resultOutPutDirectory);
		}

		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount > 0) {
			//ReverseSequenceReadToGap(scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadSet, scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount);
			GetMafftLeftAndRightSequenceTarget(scaffoldSetHead, scaffoldIndex, gapIndex, 2, resultOutPutDirectory);
			/*if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence != NULL) {
				ReverseSequence(scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence);
			}
			ReverseSequenceReadToGap(scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadSet, scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount);*/
		}
		return true;
	}
	
	return true;
}


//先获取比对的的读数
bool GetConsensusSequence(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type,char* resultOutPutDirectory) {
	char* tempSequenceFile = (char*)malloc(sizeof(char) * 450);
	strcpy(tempSequenceFile, resultOutPutDirectory);
	strcat(tempSequenceFile, "/tempSequence.fa");

	char* alignmentFile = (char*)malloc(sizeof(char) * 450);
	strcpy(alignmentFile, resultOutPutDirectory);
	strcat(alignmentFile, "/alignment.fa");

	
	ReadToGap* ReadToGap = NULL;
	long int readCount = 0;
	if (type == 0) {
		ReadToGap = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanReadSet;
		readCount = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanReadCount;
	}

	if (type == 1) {
		ReadToGap = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftReadSet;
		readCount = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftReadCount;
	}

	if (type == 2) {
		ReadToGap = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadSet;
		readCount = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount;
	}

	long int index = -1;
	long int secondIndex = -1;
	long int maxWeigth = 0;
	long int maxOverlap = 0;

	for (long int i = 0; i < readCount; i++) {
		cout << i << "----" << (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2 << "-------" << ReadToGap[i].weigth << endl;
		if (maxWeigth < ReadToGap[i].weigth) {
			index = i;
			maxWeigth = ReadToGap[i].weigth;
			maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
		}
		if (maxWeigth == ReadToGap[i].weigth && maxOverlap < (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2) {
			index = i;
			maxWeigth = ReadToGap[i].weigth;
			maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
		}
	}


	long int flag[readCount];
	for (long int i = 0; i < readCount; i++) {
		flag[i] = 1;
	}

	FILE* fp;
	if ((fp = fopen(tempSequenceFile, "w")) == NULL) {
		printf("%s, does not exist!", tempSequenceFile);
		exit(0);
	}

	for (long int i = 0; i < readCount / 2; i++) {
		secondIndex = -1;
		maxWeigth = 0;
		maxOverlap = 0;
		for (long int i = 0; i < readCount; i++) {
			if (maxWeigth < ReadToGap[i].weigth && flag[i] != 0) {
				secondIndex = i;
				maxWeigth = ReadToGap[i].weigth;
				maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
			}
			if (maxWeigth == ReadToGap[i].weigth && maxOverlap < (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2 && flag[i] != 0) {
				secondIndex = i;
				maxWeigth = ReadToGap[i].weigth;
				maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
			}

		}
		cout << "secondIndex=" << secondIndex << endl;
		fprintf(fp, ">%ld\n%s\n", secondIndex, ReadToGap[secondIndex].read);
		flag[secondIndex] = 0;
	}

	fflush(fp);
	fclose(fp);

	MAFFT(tempSequenceFile, alignmentFile);

	ContigSetHead* contigSetHead = GetContigSetFromContigSetFile(alignmentFile);

	char* consensus = GetConsensusSequenceFromOne(contigSetHead);

	int consensusLength = strlen(consensus);
	if (consensusLength > 0) {
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = (char*)malloc(sizeof(char) * (consensusLength + 1));
		strcpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence, consensus);
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence[consensusLength] = '\0';
		cout << "Mafft-Span:" << scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence << endl;
	}
	else {
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = ReadToGap[index].read;
		cout << "Mafft-Span02:" << scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence << endl;
	}

	return true;
}



bool GetMafftLeftAndRightSequenceTarget(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory) {

	ReadToGap* ReadToGap = NULL;
	long int readCount = 0;

	long int gapLength = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].gapLength;
	if (type == 0) {
		ReadToGap = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanReadSet;
		readCount = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanReadCount;
	}

	if (type == 1) {
		ReadToGap = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftReadSet;
		readCount = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftReadCount;
	}

	if (type == 2) {
		ReadToGap = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadSet;
		readCount = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount;
	}

	cout << readCount << endl;

	long int index = -1;
	long int maxOverlap = 0;
	long int maxWeigth = 0;
	long int maxReadLength = 0;
	long int diffLength = 0;

	for (long int i = 0; i < readCount; i++) {

		if (maxWeigth < ReadToGap[i].weigth) {
			index = i;
			maxWeigth = ReadToGap[i].weigth;
			maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
			maxReadLength = ReadToGap[i].readLength;

		}
		if (maxWeigth == ReadToGap[i].weigth && maxOverlap < (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2) {
			index = i;
			maxWeigth = ReadToGap[i].weigth;
			maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
			maxReadLength = ReadToGap[i].readLength;
		}
		/*if (maxWeigth == ReadToGap[i].weigth && maxOverlap == (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2 && maxReadLength < ReadToGap[i].readLength) {
			index = i;
			maxWeigth = ReadToGap[i].weigth;
			maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
			maxReadLength = ReadToGap[i].readLength;
		}*/

	}

	if (type == 0) {
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = ReadToGap[index].read;
	}
	if (type == 1) {
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence = ReadToGap[index].read;
	}
	if (type == 2) {
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence = ReadToGap[index].read;
	}

	return true;
}










//将这些读数进行比对
void MAFFT(char* tempSequenceFile, char* alignmentFile) {
	char* command = (char*)malloc(2000);
	sprintf(command, "mafft --auto --quiet --ep 1 %s > %s", tempSequenceFile, alignmentFile);
	//sprintf(command, "muscle -align %s -output %s", tempSequenceFile, alignmentFile);
	system(command);
}

//提取比对的读数信息(已有)

//根据比对信息获取一致序列
char* GetConsensusSequenceFromOne(ContigSetHead* contigSetHead) {
	
	AlignmentArrayHead* alignmentArrayHead = (AlignmentArrayHead*)malloc(sizeof(AlignmentArrayHead));
	alignmentArrayHead->alignmentArray = (AlignmentArray*)malloc(sizeof(AlignmentArray) * contigSetHead->contigSet[0].contigLength);
	for (long int i = 0; i < contigSetHead->contigSet[0].contigLength; i++) {

		alignmentArrayHead->alignmentArray[i].weight = (int*)malloc(sizeof(int) * 5);
		for (long int j = 0; j < 5; j++) {
			alignmentArrayHead->alignmentArray[i].weight[j] = 0;
		}
	}
	
	//cout << "ContigCount=" << contigSetHead->contigCount << "-contigLength" << contigSetHead->contigSet[0].contigLength << endl;
	for (long int i = 0; i < contigSetHead->contigCount; i++) {
		for (long int j = 0; j < contigSetHead->contigSet[0].contigLength; j++) {
			if (contigSetHead->contigSet[i].contig[j] == 'a' || contigSetHead->contigSet[i].contig[j] == 'A') {
				alignmentArrayHead->alignmentArray[j].weight[0]++;
			}
			if (contigSetHead->contigSet[i].contig[j] == 't' || contigSetHead->contigSet[i].contig[j] == 'T') {
				alignmentArrayHead->alignmentArray[j].weight[1]++;
			}
			if (contigSetHead->contigSet[i].contig[j] == 'g' || contigSetHead->contigSet[i].contig[j] == 'G') {
				alignmentArrayHead->alignmentArray[j].weight[2]++;
			}
			if (contigSetHead->contigSet[i].contig[j] == 'c' || contigSetHead->contigSet[i].contig[j] == 'C') {
				alignmentArrayHead->alignmentArray[j].weight[3]++;
			}
			if (contigSetHead->contigSet[i].contig[j] == '-') {
				alignmentArrayHead->alignmentArray[j].weight[4]++;
			}
		}
	}
	
	long int index = -1;
	long int position = 0;
	char* consensus = (char *)malloc(sizeof(char)* (contigSetHead->contigSet[0].contigLength+1));//此处确定加加1

	/*for (long int i = 0; i < contigSetHead->contigSet[0].contigLength; i++) {
		for (long int j = 0; j < 5; j++) {
			cout << alignmentArrayHead->alignmentArray[i].weight[j] << endl;
		}
		cout << endl;
	}*/
	
	for (long int i = 0; i < contigSetHead->contigSet[0].contigLength; i++) {
		index = VoteMajority(alignmentArrayHead->alignmentArray[i].weight, 5);
		if (index != -1) {
			if (index == 0) {
				consensus[position] = 'A';
			}
			else if (index == 1) {
				consensus[position] = 'T';
			}
			else if (index == 2) {
				consensus[position] = 'G';
			}
			else {
				consensus[position] = 'C';
			}
		}
		if (index != -1) {
			position++;
		}

	}
	
	consensus[position] = '\0';
	return consensus;

}

long int VoteMajority(int* weight, int count) {
	int max = 0;
	long int index = -1;
	for (long int i = 0; i < count; i++) {
		if (weight[i] > max) {
			max = weight[i];
			index = i;
		}
	}

	max = 0;
	long int index1 = -1;
	for (long int i = 0; i < count; i++) {
		if (i == index) {
			continue;
		}
		if (weight[i] > max) {
			max = weight[i];
			index1 = i;
		}
	}

	if (index != count - 1) {//index！=4
		return index;
		if ((double)weight[index1] / weight[index] < 0.6) {
			return index;
		}
	}

	return -1;


}

void ReverseSequenceReadToGap(ReadToGap* readToGap, long int count) {
	for (long int p = 0; p < count; p++) {
		ReverseSequence(readToGap[p].read);
	}
}




#endif