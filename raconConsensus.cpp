#ifndef RACONCONSENSUS_CPP_INCLUDED 
#define RACONCONSENSUS_CPP_INCLUDED 

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

#include "raconConsensus.h"


using namespace std;

bool GetRaconConsensusSequenceForOneGap(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, char* resultOutPutDirectory) {
	

	if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanReadCount > 0) {
		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].estimatedGapDistance > 100) {
			GetSpanSequenceTarget(scaffoldSetHead, scaffoldIndex, gapIndex,0, resultOutPutDirectory);
		}
		if ((scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].estimatedGapDistance < 100 && scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].estimatedGapDistance>0) || scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence == NULL) {
			GetConsenSequence(scaffoldSetHead, scaffoldIndex, gapIndex, 0, resultOutPutDirectory);
		}
		
		return true;
	}
	else {
		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftReadCount > 0) {
			GetleftAndRightConsenSequence(scaffoldSetHead, scaffoldIndex, gapIndex, 1, resultOutPutDirectory);

			/*if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence == NULL) {
				GetConsenSequence(scaffoldSetHead, scaffoldIndex, gapIndex, 1, resultOutPutDirectory);
			}*/
		}
		
		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount > 0) {

			ReverseSequenceReadToGap(scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadSet, scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount);

			GetleftAndRightConsenSequence(scaffoldSetHead, scaffoldIndex, gapIndex, 2, resultOutPutDirectory);
			
			if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence != NULL) {
				ReverseSequence(scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence);
			}
			
			ReverseSequenceReadToGap(scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadSet, scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightReadCount);

			/*if (scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence == NULL) {
				GetConsenSequence(scaffoldSetHead, scaffoldIndex, gapIndex, 2, resultOutPutDirectory);
			}*/
		}
		return true;
	}
	
	//cout << "GetRaconConsensusSequenceForOneGap" << endl;
	return false;
}




//针对左右端读数的一致序列
bool GetleftAndRightConsenSequence(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory) {
	char* targetFile = (char*)malloc(sizeof(char) * 450);
	strcpy(targetFile, resultOutPutDirectory);
	strcat(targetFile, "/targetFile.fa");

	char* sequenceFile = (char*)malloc(sizeof(char) * 450);
	strcpy(sequenceFile, resultOutPutDirectory);
	strcat(sequenceFile, "/sequenceFile.fa");

	char* overlapFile = (char*)malloc(sizeof(char) * 450);
	strcpy(overlapFile, resultOutPutDirectory);
	strcat(overlapFile, "/overlapFile.paf");

	char* consensusFile = (char*)malloc(sizeof(char) * 450);
	strcpy(consensusFile, resultOutPutDirectory);
	strcat(consensusFile, "/consensusFile.fa");

	char* dataBase = (char*)malloc(sizeof(char) * 450);
	strcpy(dataBase, resultOutPutDirectory);
	strcat(dataBase, "/dataBase.fa");

	char* blastnDataBase = (char*)malloc(sizeof(char) * 450);
	strcpy(blastnDataBase, resultOutPutDirectory);
	strcat(blastnDataBase, "/blastn/database");

	char* blastnResult = (char*)malloc(sizeof(char) * 450);
	strcpy(blastnResult, resultOutPutDirectory);
	strcat(blastnResult, "/blastn/result");


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

	ReadSort* readSort = NULL;
	if (NULL == (readSort = (ReadSort*)malloc(sizeof(ReadSort)))) {
		perror("readSort malloc error");
		exit(1);
	}
	readSort->readSortForLength = NULL;
	readSort->sortLengthCount = 0;
	readSort->readSortForLength = (ReadSortForLength*)malloc(sizeof(ReadSortForLength) * 20);

	for (long int i = 0; i < 20; i++) {//分为20类
		readSort->readSortForLength[i].readSortCountForLengthSet = (ReadSortCountForLengthSet*)malloc(sizeof(ReadSortCountForLengthSet) * 500);
		readSort->readSortForLength[i].readSortForLengthCount = 0;
		readSort->readSortForLength[i].aveReadLength = 0;

		readSort->readSortForLength[i].readSortCountForSimilar = (ReadSortForSimilar*)malloc(sizeof(ReadSortForSimilar) * 20);
		readSort->readSortForLength[i].sortSimilarCount = 0;
		for (long int j = 0; j < 500; j++) {//每类200个
			readSort->readSortForLength[i].readSortCountForLengthSet[j].readIndex = -1;
			readSort->readSortForLength[i].readSortCountForLengthSet[j].flag = 0;
		}

		for (long int k = 0; k < 20; k++) {//
			readSort->readSortForLength[i].readSortCountForSimilar[k].readSortCountForSimilarSet = (ReadSortCountForSimilarSet*)malloc(sizeof(ReadSortCountForSimilarSet) * 500);
			readSort->readSortForLength[i].readSortCountForSimilar[k].readSortForSimilarCount = 0;
			for (long int m = 0; m < 500; m++) {
				readSort->readSortForLength[i].readSortCountForSimilar[k].readSortCountForSimilarSet[m].readIndex = -1;
				readSort->readSortForLength[i].readSortCountForSimilar[k].readSortCountForSimilarSet[m].flag = 0;

			}

		}
	}

	long int index = -1;
	long int maxOverlap = 0;
	long int maxWeigth = 0;
	long int maxReadLength = 0;
	long int readSortCount = 1;
	long int flagIndex = 0;
	long int sortNum = 0;
	long int sortReadNum = 0;

	//
	for (long int i = 0; i < readCount; i++) {

		//cout << "i=" << i << "   readLength=" << ReadToGap[i].readLength << "  flagLength=" << ReadToGap[flagIndex].readLength << "  diff=" << 0.3 * ReadToGap[flagIndex].readLength << endl;
		if (abs(ReadToGap[flagIndex].readLength - ReadToGap[i].readLength) > 0.5 * ReadToGap[flagIndex].readLength) {
			readSortCount = readSortCount + 1;
			flagIndex = i;
		}
	}
	//cout << endl;
	//cout << "readSortCount=" << readSortCount << endl;
	
	//按长度分类
	flagIndex = 0;
	for (long int i = 0; i < readCount; i++) {
		if (abs(ReadToGap[flagIndex].readLength - ReadToGap[i].readLength) > 0.5 * ReadToGap[flagIndex].readLength) {
			readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength / readSort->readSortForLength[sortNum].readSortForLengthCount;

			sortNum = sortNum + 1;
			sortReadNum = 0;
			readSort->readSortForLength[sortNum].readSortCountForLengthSet[sortReadNum].readIndex = i;
			readSort->readSortForLength[sortNum].readSortForLengthCount = readSort->readSortForLength[sortNum].readSortForLengthCount + 1;
			readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength + ReadToGap[i].readLength;
			flagIndex = i;
			sortReadNum = sortReadNum + 1;
		}
		else {
			//cout << "i=" << i << "   readLength=" << ReadToGap[i].readLength << "  flagLength=" << ReadToGap[flagIndex].readLength << "  diff=" << 0.3 * ReadToGap[flagIndex].readLength << endl;

			readSort->readSortForLength[sortNum].readSortCountForLengthSet[sortReadNum].readIndex = i;

			readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength + ReadToGap[i].readLength;
			readSort->readSortForLength[sortNum].readSortForLengthCount = readSort->readSortForLength[sortNum].readSortForLengthCount + 1;
			sortReadNum = sortReadNum + 1;

		}
	}
	readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength / readSort->readSortForLength[sortNum].readSortForLengthCount;
	readSort->sortLengthCount = sortNum + 1;


	/*
	cout << "sortLengthCount=" << readSort->sortLengthCount << endl;
	for (long int i = 0; i < readSort->sortLengthCount; i++) {
		cout << "aveReadLength=" << readSort->readSortForLength[i].aveReadLength << "    readSortForLengthCount=" << readSort->readSortForLength[i].readSortForLengthCount << endl;
		for (long int j = 0; j < readSort->readSortForLength[i].readSortForLengthCount; j++) {
			cout << "readIndex=" << readSort->readSortForLength[i].readSortCountForLengthSet[j].readIndex << endl;
		}
	}*/
	FILE* fp;
	char command[3000];
	//按相似性分类
	for (long int i = 0; i < readSort->sortLengthCount; i++) {
		if ((fp = fopen(dataBase, "w")) == NULL) {
			printf("%s, does not exist!", dataBase);
			exit(0);
		}
		for (long int j = 0; j < readSort->readSortForLength[i].readSortForLengthCount; j++) {//
			fprintf(fp, ">%ld\n%s\n", j, ReadToGap[readSort->readSortForLength[i].readSortCountForLengthSet[j].readIndex].read);
		}
		fflush(fp);
		fclose(fp);

		sprintf(command, "makeblastdb -in %s -dbtype nucl -out %s", dataBase, blastnDataBase);
		system(command);
		sprintf(command, "blastn -query %s -db %s -out %s -evalue 1e-10 -outfmt 6 -perc_identity 80 -max_hsps 1", dataBase, blastnDataBase, blastnResult);//-max_hsps 1
		system(command);
		const char* split = "\t";
		char* p;
		long int qIndex;
		long int tIndex;
		long int startIndex = -1;
		long int flag = 0;
		long int maxSize = 1000;
		long int num[readSort->readSortForLength[i].readSortForLengthCount];
		for (long int j = 0; j < readSort->readSortForLength[i].readSortForLengthCount; j++) {
			num[j] = 1;
		}
		long int sortForSimilarCount = 0;

		char* line = (char*)malloc(sizeof(char) * maxSize);

		if ((fp = fopen(blastnResult, "r")) == NULL) {
			printf("%s, does not exist!", blastnResult);
			exit(0);
		}
		sortNum = -1;
		sortReadNum = 0;

		while ((fgets(line, maxSize, fp)) != NULL) {//
			p = strtok(line, split);
			qIndex = atoi(p);
			//cout << "qname=" << qIndex << "----tname=";
			p = strtok(NULL, split);
			tIndex = atoi(p);
			//cout << tIndex << "    startIndex="<< startIndex<<endl;

			if (startIndex != qIndex && num[qIndex] == 0) {
				//cout << "butianjia" << endl;
				flag = 0;
			}

			if (startIndex != qIndex && num[qIndex] == 1) {
				//cout << "tianjia*****-*"<< endl;
				//cout << num[qIndex] << endl;
				readSort->readSortForLength[i].sortSimilarCount = readSort->readSortForLength[i].sortSimilarCount + 1;//
				sortNum = sortNum + 1;
				sortReadNum = 0;

				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortCountForSimilarSet[sortReadNum].readIndex = readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex;
				//cout << readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex << endl;
				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount = readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount + 1;//
				sortReadNum = sortReadNum + 1;
				num[tIndex] = 0;
				flag = 1;
			}


			if (startIndex == qIndex && num[tIndex] == 1 && flag == 1) {
				//cout << "tianjia******"<< endl;

				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortCountForSimilarSet[sortReadNum].readIndex = readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex;
				//cout << readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex << endl;
				sortReadNum = sortReadNum + 1;
				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount = readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount + 1;
				num[tIndex] = 0;
			}
			startIndex = qIndex;
		}


		//
		/*
		cout << "twoSortCount=" << readSort->readSortForLength[i].sortSimilarCount << endl;
		for (long int j = 0; j < readSort->readSortForLength[i].sortSimilarCount; j++) {
			cout << "******twoSortCount-readCount=" << readSort->readSortForLength[i].readSortCountForSimilar[j].readSortForSimilarCount << endl;
			for (long int k = 0; k < readSort->readSortForLength[i].readSortCountForSimilar[j].readSortForSimilarCount; k++) {
				cout << "readIndex=" << readSort->readSortForLength[i].readSortCountForSimilar[j].readSortCountForSimilarSet[k].readIndex << endl;
			}
		}*/
	}

	fflush(fp);
	fclose(fp);
	long int firstSortIndex = -1;
	long int twoSortIndex = -1;
	long int maxCount = 0;
	long int targetIndex = -1;

	//cout << "*************************************************consensus****************************************************************" << endl;
	//cout << endl;
	
	//
	//cout << "sortCount=" << readSort->sortLengthCount << endl;
	//cout << "gapLength=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].gapLength << "  estimatedGapDistance=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].estimatedGapDistance << endl;

	if (readSort->sortLengthCount == 1) {//
		firstSortIndex = 0;
		//cout << "aveReadLength=" << readSort->readSortForLength[0].aveReadLength << endl;
	}
	else {
		for (long int i = 0; i < readSort->sortLengthCount; i++) {//
			//cout << "readSortForLengthCount=" << readSort->readSortForLength[i].readSortForLengthCount << endl;
			if (readSort->readSortForLength[i].readSortForLengthCount> maxCount) {
				maxCount = readSort->readSortForLength[i].readSortForLengthCount;
				firstSortIndex = i;
			}
		}

		if (firstSortIndex == -1) {//
			firstSortIndex = 0;
		}
	}


	//cout << "firstSortIndex=" << firstSortIndex << endl;
	//cout << "twoSortCount=" << readSort->readSortForLength[firstSortIndex].sortSimilarCount << endl;
	
	if (readSort->readSortForLength[firstSortIndex].sortSimilarCount == 1) {//
		twoSortIndex = 0;
	}
	else {

		for (long int j = 0; j < readSort->readSortForLength[firstSortIndex].sortSimilarCount; j++) {
			//cout << "2leireadCount=" << readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[j].readSortForSimilarCount << endl;
			if (readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[j].readSortForSimilarCount > maxCount) {
				maxCount = readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[j].readSortForSimilarCount;
				twoSortIndex = j;
			}
		}
	}
	//cout << "twoSortIndex=" << twoSortIndex << endl;

	if (firstSortIndex != -1 && twoSortIndex != -1) {

		long int realReadIndex = -1;
		long int index = -1;
		long int maxWeigth = 0;
		long int maxOverlap = 0;
		//cout << "twoReadConut=" << readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortForSimilarCount << endl;
		for (long int k = 0; k < readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortForSimilarCount; k++) {
			realReadIndex = readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortCountForSimilarSet[k].readIndex;
			//cout << "realReadIndex=" << realReadIndex << endl;
			if (maxWeigth < ReadToGap[realReadIndex].weigth) {
				index = realReadIndex;
				maxWeigth = ReadToGap[realReadIndex].weigth;
				maxOverlap = (ReadToGap[realReadIndex].leftContigOverlap + ReadToGap[realReadIndex].rightContigOverlap) / 2;
			}
			if (maxWeigth == ReadToGap[realReadIndex].weigth && maxOverlap < (ReadToGap[realReadIndex].leftContigOverlap + ReadToGap[realReadIndex].rightContigOverlap) / 2) {
				index = realReadIndex;
				maxWeigth = ReadToGap[realReadIndex].weigth;
				maxOverlap = (ReadToGap[realReadIndex].leftContigOverlap + ReadToGap[realReadIndex].rightContigOverlap) / 2;
			}
		}
		//cout << "index=" << index << endl;
		if ((fp = fopen(targetFile, "w")) == NULL) {
			printf("%s, does not exist!", targetFile);
			exit(0);
		}
		fprintf(fp, ">%ld\n%s\n", index, ReadToGap[index].read);
		fflush(fp);
		fclose(fp);


		if ((fp = fopen(sequenceFile, "w")) == NULL) {
			printf("%s, does not exist!", sequenceFile);
			exit(0);
		}
		for (long int k = 0; k < readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortForSimilarCount; k++) {
			realReadIndex = readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortCountForSimilarSet[k].readIndex;
			if (realReadIndex != index) {
				fprintf(fp, ">%ld\n%s\n", realReadIndex, ReadToGap[realReadIndex].read);
			}

		}
		fflush(fp);
		fclose(fp);

		sprintf(command, "minimap2 -k14 -w5 -n2 -m20 -s 40 --sr --frag yes %s %s > %s", targetFile, sequenceFile, overlapFile);
		system(command);

		sprintf(command, "racon -t10 %s %s %s > %s 2>racon.log", sequenceFile, overlapFile, targetFile, consensusFile);
		system(command);


		if (type == 1) {
			if (access(consensusFile, F_OK) != 0) {
				scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence = ReadToGap[index].read;

				//cout << "**************L01********************************" << endl;
				return true;
			}
			if (FileIsNull(consensusFile) == 1) {
				scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence = ReadToGap[index].read;

				//cout << "**************L02********************************" << endl;
				return true;
			}
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence = GetConsensusSequenceFromFile(consensusFile);
			//cout << scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence << endl;
			return true;
		}

		if (type == 2) {
			if (access(consensusFile, F_OK) != 0) {
				scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence = ReadToGap[index].read;

				//cout << "**************R01********************************" << endl;
				return true;
			}
			if (FileIsNull(consensusFile) == 1) {
				scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence = ReadToGap[index].read;

				//cout << "**************R02********************************" << endl;
				return true;
			}
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence = GetConsensusSequenceFromFile(consensusFile);
			//cout << scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence << endl;
			return true;
		}
	}
}


//针对跨越读数读数的一致序列

bool GetSpanSequenceTarget(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type,char *resultOutPutDirectory) {

	char* targetFile = (char*)malloc(sizeof(char) * 450);
	strcpy(targetFile, resultOutPutDirectory);
	strcat(targetFile, "/targetFile.fa");

	char* sequenceFile = (char*)malloc(sizeof(char) * 450);
	strcpy(sequenceFile, resultOutPutDirectory);
	strcat(sequenceFile, "/sequenceFile.fa");

	char* overlapFile = (char*)malloc(sizeof(char) * 450);
	strcpy(overlapFile, resultOutPutDirectory);
	strcat(overlapFile, "/overlapFile.paf");

	char* consensusFile = (char*)malloc(sizeof(char) * 450);
	strcpy(consensusFile, resultOutPutDirectory);
	strcat(consensusFile, "/consensusFile.fa");

	char* dataBase = (char*)malloc(sizeof(char) * 450);
	strcpy(dataBase, resultOutPutDirectory);
	strcat(dataBase, "/dataBase.fa");

	char* blastnDataBase = (char*)malloc(sizeof(char) * 450);
	strcpy(blastnDataBase, resultOutPutDirectory);
	strcat(blastnDataBase, "/blastn/database");

	char* blastnResult = (char*)malloc(sizeof(char) * 450);
	strcpy(blastnResult, resultOutPutDirectory);
	strcat(blastnResult, "/blastn/result");
	
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

	//cout << readCount << endl;

	ReadSort* readSort = NULL;
	if (NULL == (readSort = (ReadSort*)malloc(sizeof(ReadSort)))) {
		perror("readSort malloc error");
		exit(1);
	}
	readSort->readSortForLength = NULL;
	readSort->sortLengthCount = 0;
	readSort->readSortForLength = (ReadSortForLength*)malloc(sizeof(ReadSortForLength) * 20);
	for (long int i = 0; i < 20; i++) {
		readSort->readSortForLength[i].readSortCountForLengthSet = (ReadSortCountForLengthSet*)malloc(sizeof(ReadSortCountForLengthSet) * 500);
		readSort->readSortForLength[i].readSortForLengthCount = 0;
		readSort->readSortForLength[i].aveReadLength = 0;

		readSort->readSortForLength[i].readSortCountForSimilar = (ReadSortForSimilar*)malloc(sizeof(ReadSortForSimilar) * 20);
		readSort->readSortForLength[i].sortSimilarCount = 0;
		for (long int j = 0; j < 500; j++) {
			readSort->readSortForLength[i].readSortCountForLengthSet[j].readIndex = -1;
			readSort->readSortForLength[i].readSortCountForLengthSet[j].flag = 0;
			
		}

		for (long int k = 0; k < 20; k++) {
			readSort->readSortForLength[i].readSortCountForSimilar[k].readSortCountForSimilarSet = (ReadSortCountForSimilarSet*)malloc(sizeof(ReadSortCountForSimilarSet) * 500);
			readSort->readSortForLength[i].readSortCountForSimilar[k].readSortForSimilarCount = 0;
			for (long int m = 0; m < 500; m++) {
				readSort->readSortForLength[i].readSortCountForSimilar[k].readSortCountForSimilarSet[m].readIndex = -1;
				readSort->readSortForLength[i].readSortCountForSimilar[k].readSortCountForSimilarSet[m].flag = 0;

			}

		}
	}

	long int index = -1;
	long int maxOverlap = 0;
	long int maxWeigth = 0;
	long int maxReadLength = 0;
	long int readSortCount = 1;
	long int flagIndex = 0;
	long int sortNum = 0;
	long int sortReadNum = 0;

	//
	for (long int i = 0; i < readCount; i++) {

		//cout << "i=" << i << "   readLength=" << ReadToGap[i].readLength << "  flagLength=" << ReadToGap[flagIndex].readLength <<"  diff="<< 0.3 * ReadToGap[flagIndex].readLength << endl;
		if (abs(ReadToGap[flagIndex].readLength - ReadToGap[i].readLength) > 0.3 * ReadToGap[flagIndex].readLength) {
			readSortCount = readSortCount + 1;
			flagIndex = i;
		}
	}
	//cout << endl;
	//cout << "readSortCount=" << readSortCount << endl;

	//
	flagIndex = 0;
	for (long int i = 0; i < readCount; i++) {
		if (abs(ReadToGap[flagIndex].readLength - ReadToGap[i].readLength) > 0.3 * ReadToGap[flagIndex].readLength) {
			readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength / readSort->readSortForLength[sortNum].readSortForLengthCount;

			sortNum = sortNum + 1;
			sortReadNum = 0;
			readSort->readSortForLength[sortNum].readSortCountForLengthSet[sortReadNum].readIndex =i;
			readSort->readSortForLength[sortNum].readSortForLengthCount = readSort->readSortForLength[sortNum].readSortForLengthCount + 1;
			readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength + ReadToGap[i].readLength;
			flagIndex = i;
			sortReadNum = sortReadNum + 1;
		}
		else {
			//cout << "i=" << i << "   readLength=" << ReadToGap[i].readLength << "  flagLength=" << ReadToGap[flagIndex].readLength << "  diff=" << 0.3 * ReadToGap[flagIndex].readLength << endl;
			
			readSort->readSortForLength[sortNum].readSortCountForLengthSet[sortReadNum].readIndex = i;

			readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength + ReadToGap[i].readLength;
			readSort->readSortForLength[sortNum].readSortForLengthCount = readSort->readSortForLength[sortNum].readSortForLengthCount + 1;
			sortReadNum = sortReadNum + 1;
			
		}
	}
	readSort->readSortForLength[sortNum].aveReadLength = readSort->readSortForLength[sortNum].aveReadLength / readSort->readSortForLength[sortNum].readSortForLengthCount;
	readSort->sortLengthCount = sortNum+1;

	
	//
	//cout << "sortLengthCount=" << readSort->sortLengthCount << endl;
	for (long int i = 0; i < readSort->sortLengthCount; i++) {
		//cout << "aveReadLength=" << readSort->readSortForLength[i].aveReadLength << "    readSortForLengthCount=" << readSort->readSortForLength[i].readSortForLengthCount << endl;
		for (long int j = 0; j < readSort->readSortForLength[i].readSortForLengthCount; j++) {
			//cout << "readIndex=" << readSort->readSortForLength[i].readSortCountForLengthSet[j].readIndex << endl;
		}
	}
	FILE* fp;
	char command[3000];
	//
	for (long int i = 0; i < readSort->sortLengthCount; i++) {//
		if ((fp = fopen(dataBase, "w")) == NULL) {
			printf("%s, does not exist!", dataBase);
			exit(0);
		}
		for (long int j = 0; j < readSort->readSortForLength[i].readSortForLengthCount; j++) {//
			fprintf(fp, ">%ld\n%s\n", j, ReadToGap[readSort->readSortForLength[i].readSortCountForLengthSet[j].readIndex].read);
		}
		fflush(fp);
		fclose(fp);

		sprintf(command, "makeblastdb -in %s -dbtype nucl -out %s", dataBase, blastnDataBase);
		system(command);
		sprintf(command, "blastn -query %s -db %s -out %s -evalue 1e-10 -outfmt 6 -perc_identity 80", dataBase, blastnDataBase, blastnResult);//-max_hsps 1
		system(command);
		const char* split = "\t";
		char* p;
		long int qIndex;
		long int tIndex;
		long int startIndex=-1;
		long int flag = 0;
		long int maxSize = 1000;
		long int num[readSort->readSortForLength[i].readSortForLengthCount];
		for (long int j = 0; j < readSort->readSortForLength[i].readSortForLengthCount; j++) {
			num[j] = 1;
		}
		long int sortForSimilarCount = 0;

		char* line = (char*)malloc(sizeof(char) * maxSize);

		if ((fp = fopen(blastnResult, "r")) == NULL) {
			printf("%s, does not exist!", blastnResult);
			exit(0);
		}
		sortNum = -1;
		sortReadNum = 0;
		
		while ((fgets(line, maxSize, fp)) != NULL) {//
			p = strtok(line, split);
			qIndex = atoi(p);
			//cout << "qname=" << qIndex << "----tname=";
			p = strtok(NULL, split);
			tIndex = atoi(p);
			//cout << tIndex << "    startIndex="<< startIndex<<endl;

			if (startIndex != qIndex && num[qIndex] == 0) {
				//cout << "butianjia" << endl;
				flag = 0;
			}
			
			if (startIndex!= qIndex && num[qIndex] == 1) {
				//cout << "tianjia*****-*"<< endl;
				//cout << num[qIndex] << endl;
				readSort->readSortForLength[i].sortSimilarCount = readSort->readSortForLength[i].sortSimilarCount + 1;
				sortNum = sortNum + 1;
				sortReadNum = 0;

				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortCountForSimilarSet[sortReadNum].readIndex = readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex;
				//cout << readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex << endl;
				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount = readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount + 1;
				sortReadNum = sortReadNum + 1;
				num[tIndex] = 0;
				flag = 1;
			}
			

			if(startIndex == qIndex && num[tIndex] == 1 && flag==1){
				//cout << "tianjia******"<< endl;
				
				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortCountForSimilarSet[sortReadNum].readIndex = readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex;
				//cout << readSort->readSortForLength[i].readSortCountForLengthSet[tIndex].readIndex << endl;
				sortReadNum = sortReadNum + 1;
				readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount = readSort->readSortForLength[i].readSortCountForSimilar[sortNum].readSortForSimilarCount + 1;
				num[tIndex] = 0;
			}
			startIndex = qIndex;
		}
		

        //
		//cout << "twoSortCount=" << readSort->readSortForLength[i].sortSimilarCount << endl;
		for (long int j = 0; j < readSort->readSortForLength[i].sortSimilarCount; j++) {
			//cout << "******twoSortCount-readCount=" << readSort->readSortForLength[i].readSortCountForSimilar[j].readSortForSimilarCount << endl;
			for (long int k = 0; k < readSort->readSortForLength[i].readSortCountForSimilar[j].readSortForSimilarCount; k++) {
				//cout << "readIndex=" << readSort->readSortForLength[i].readSortCountForSimilar[j].readSortCountForSimilarSet[k].readIndex << endl;
			}
		}
	}

	fflush(fp);
	fclose(fp);
	long int firstSortIndex = -1;
	long int twoSortIndex = -1;
	long int maxCount = 0;
	long int targetIndex = -1;

	//cout << "*************************************************consensus****************************************************************" << endl;
	//cout << endl;

	
	//
	//cout << "sortCount=" << readSort->sortLengthCount << endl;
	//cout << "gapLength=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].gapLength << "  estimatedGapDistance=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].estimatedGapDistance << endl;
	
	if (readSort->sortLengthCount==1) {//
		firstSortIndex = 0;
		//cout << "aveReadLength=" << readSort->readSortForLength[0].aveReadLength << endl;
	}
	else {

		for (long int i = 0; i < readSort->sortLengthCount; i++) {//
			//cout << "readSortForLengthCount=" << readSort->readSortForLength[i].readSortForLengthCount << endl;
			if (readSort->readSortForLength[i].readSortForLengthCount > maxCount) {
				maxCount = readSort->readSortForLength[i].readSortForLengthCount;
				firstSortIndex = i;
			}
		}

		//for (long int i = 0; i < readSort->sortLengthCount; i++) {//
		//	cout << "aveReadLength=" << readSort->readSortForLength[i].aveReadLength << endl;
		//	if (abs(readSort->readSortForLength[i].aveReadLength - scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].gapLength) < 0.3 * scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].gapLength) {
		//		firstSortIndex = i;
		//	}
		//}

		//if (firstSortIndex == -1) {//
		//	for (long int i = 0; i < readSort->sortLengthCount; i++) {
		//		if (abs(readSort->readSortForLength[i].aveReadLength - scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].estimatedGapDistance) < 0.3 * scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].estimatedGapDistance) {
		//			firstSortIndex = i;
		//		}
		//	}
		//}

		if (firstSortIndex == -1) {//
			firstSortIndex = 0;
		}
	}
	
	maxCount = 0;
	//cout << "firstSortIndex=" << firstSortIndex << endl;
	//cout << "twoSortCount=" << readSort->readSortForLength[firstSortIndex].sortSimilarCount << endl;
	if (readSort->readSortForLength[firstSortIndex].sortSimilarCount == 1) {//
		twoSortIndex = 0;
	}
	else {
		
		for (long int j = 0; j < readSort->readSortForLength[firstSortIndex].sortSimilarCount; j++) {
			//cout << "2leireadCount=" << readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[j].readSortForSimilarCount << endl;
			if (readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[j].readSortForSimilarCount > maxCount) {
				maxCount = readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[j].readSortForSimilarCount;
				twoSortIndex = j;
			}
		}
	}
	//cout << "twoSortIndex=" << twoSortIndex << endl;

	if (firstSortIndex != -1 && twoSortIndex!=-1) {
		
		long int realReadIndex = -1;
		long int index = -1;
		long int maxWeigth = 0;
		long int maxOverlap = 0;
		//cout << "twoReadConut=" << readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortForSimilarCount << endl;
		for (long int k = 0; k < readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortForSimilarCount; k++) {
			realReadIndex = readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortCountForSimilarSet[k].readIndex;
			//cout << "realReadIndex=" << realReadIndex << endl;
			if (maxWeigth < ReadToGap[realReadIndex].weigth) {
				index = realReadIndex;
				maxWeigth = ReadToGap[realReadIndex].weigth;
				maxOverlap = (ReadToGap[realReadIndex].leftContigOverlap + ReadToGap[realReadIndex].rightContigOverlap) / 2;
			}
			if (maxWeigth == ReadToGap[realReadIndex].weigth && maxOverlap < (ReadToGap[realReadIndex].leftContigOverlap + ReadToGap[realReadIndex].rightContigOverlap) / 2) {
				index = realReadIndex;
				maxWeigth = ReadToGap[realReadIndex].weigth;
				maxOverlap = (ReadToGap[realReadIndex].leftContigOverlap + ReadToGap[realReadIndex].rightContigOverlap) / 2;
			}
		}
		//cout << "index=" << index << endl;
		if ((fp = fopen(targetFile, "w")) == NULL) {
			printf("%s, does not exist!", targetFile);
			exit(0);
		}
		fprintf(fp, ">%ld\n%s\n", index, ReadToGap[index].read);
		fflush(fp);
		fclose(fp);


		if ((fp = fopen(sequenceFile, "w")) == NULL) {
			printf("%s, does not exist!", sequenceFile);
			exit(0);
		}
		for (long int k = 0; k < readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortForSimilarCount; k++) {
			realReadIndex = readSort->readSortForLength[firstSortIndex].readSortCountForSimilar[twoSortIndex].readSortCountForSimilarSet[k].readIndex;
			if (realReadIndex != index) {
				fprintf(fp, ">%ld\n%s\n", realReadIndex, ReadToGap[realReadIndex].read);
			}

		}
		fflush(fp);
		fclose(fp);

		sprintf(command, "minimap2 -k14 -w5 -n2 -m20 -s 40 --sr --frag yes %s %s > %s", targetFile, sequenceFile, overlapFile);
		system(command);

		sprintf(command, "racon -t10 %s %s %s > %s 2>racon.log", sequenceFile, overlapFile, targetFile, consensusFile);
		system(command);

		if (type == 0) {
			if (access(consensusFile, F_OK) != 0) {
				scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = ReadToGap[index].read;

				//cout << "**************01********************************" << endl;
				return true;
			}
			if (FileIsNull(consensusFile) == 1) {
				scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = ReadToGap[index].read;

				//cout << "**************02********************************" << endl;
				return true;
			}
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = GetConsensusSequenceFromFile(consensusFile);
			//cout << scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence << endl;
			return true;
		}
	}
	return false;
}

//处理GAP比较短的的gap
bool GetConsenSequence(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory) {

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
		//if (maxWeigth == ReadToGap[i].weigth && maxOverlap == (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2 && maxReadLength < ReadToGap[i].readLength) {
		//	index = i;
		//	maxWeigth = ReadToGap[i].weigth;
		//	maxOverlap = (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2;
		//	maxReadLength = ReadToGap[i].readLength;
		//}

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


//ͨ
/*
bool GetSpanSequenceTarget(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex, long int gapIndex, int type, char* resultOutPutDirectory) {

	char* targetFile = (char*)malloc(sizeof(char) * 450);
	strcpy(targetFile, resultOutPutDirectory);
	strcat(targetFile, "/targetFile.fa");

	char* sequenceFile = (char*)malloc(sizeof(char) * 450);
	strcpy(sequenceFile, resultOutPutDirectory);
	strcat(sequenceFile, "/sequenceFile.fa");

	char* overlapFile = (char*)malloc(sizeof(char) * 450);
	strcpy(overlapFile, resultOutPutDirectory);
	strcat(overlapFile, "/overlapFile.paf");

	char* consensusFile = (char*)malloc(sizeof(char) * 450);
	strcpy(consensusFile, resultOutPutDirectory);
	strcat(consensusFile, "/consensusFile.fa");

	char* consensusFile1 = (char*)malloc(sizeof(char) * 450);
	strcpy(consensusFile1, resultOutPutDirectory);
	strcat(consensusFile1, "/consensusFile1.fa");

	char* dataBase = (char*)malloc(sizeof(char) * 450);
	strcpy(dataBase, resultOutPutDirectory);
	strcat(dataBase, "/dataBase.fa");

	char* blastnDataBase = (char*)malloc(sizeof(char) * 450);
	strcpy(blastnDataBase, resultOutPutDirectory);
	strcat(blastnDataBase, "/blastn/database");

	char* blastnResult = (char*)malloc(sizeof(char) * 450);
	strcpy(blastnResult, resultOutPutDirectory);
	strcat(blastnResult, "/blastn/result");

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

	cout << readCount << endl;

	long int index = -1;
	long int maxOverlap = 0;
	long int maxWeigth = 0;
	long int maxReadLength = 0;

	FILE* fp;
	if ((fp = fopen(dataBase, "w")) == NULL) {
		printf("%s, does not exist!", dataBase);
		exit(0);
	}

	for (long int i = 0; i < readCount; i++) {
		//cout << i << "----" << (ReadToGap[i].leftContigOverlap + ReadToGap[i].rightContigOverlap) / 2 << "-------" << ReadToGap[i].weigth << endl;
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
		fprintf(fp, ">%ld\n%s\n", i, ReadToGap[i].read);
	}
	cout << "scaffoldIndex=" << scaffoldIndex << "   gapIndex=" << gapIndex << endl;
	cout << "index=" << index <<"---diff:"<< ReadToGap[index].leftContigDiff<<"-----"<< ReadToGap[index].rightContigDiff<< endl;
	fflush(fp);
	fclose(fp);

	/*
	if ((fp = fopen(targetFile, "w")) == NULL) {
		printf("%s, does not exist!", targetFile);
		exit(0);
	}
	fprintf(fp, ">%ld\n%s\n", index, ReadToGap[index].read);
	fflush(fp);
	fclose(fp);

	if ((fp = fopen(sequenceFile, "w")) == NULL) {
		printf("%s, does not exist!", sequenceFile);
		exit(0);
	}
	long int flag[readCount];
	for (long int i = 0; i < readCount; i++) {
		flag[i] = 1;
	}
	for (long int j = 0; j < 5; j++) {
		for (long int i = 0; i < readCount; i++) {
			if (abs(ReadToGap[i].readLength - ReadToGap[index].readLength) < 1000 && index != i) {
				fprintf(fp, ">%ld\n%s\n", i, ReadToGap[i].read);
				flag[i] = 0;

			}
		}
	}


	fflush(fp);
	fclose(fp);

	char command[3000];


	sprintf(command, "minimap2 -k14 -w5 -n2 -m20 -s 40 --sr --frag yes %s %s > %s", targetFile, sequenceFile, overlapFile);
	system(command);

	sprintf(command, "racon -t10 %s %s %s > %s 2>racon.log", sequenceFile, overlapFile, targetFile, consensusFile);
	system(command);

	if (type == 0) {
		if (access(consensusFile, F_OK) != 0) {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = ReadToGap[index].read;

			cout << "**************01********************************" << endl;
			return true;
		}
		if (FileIsNull(consensusFile) == 1) {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = ReadToGap[index].read;

			cout << "**************01********************************" << endl;
			return true;
		}
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence = GetConsensusSequenceFromFile(consensusFile);
		cout << scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].spanConsensusSequence << endl;


	}
	if (type == 1) {
		if (access(consensusFile, F_OK) != 0) {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence = ReadToGap[index].read;
			return true;
		}
		if (FileIsNull(consensusFile) == 1) {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence = ReadToGap[index].read;
			return true;
		}
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence = GetConsensusSequenceFromFile(consensusFile);
		//cout << scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].leftConsensusSequence << endl;
	}
	if (type == 2) {
		if (access(consensusFile, F_OK) != 0) {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence = ReadToGap[index].read;
			return true;
		}
		if (FileIsNull(consensusFile) == 1) {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence = ReadToGap[index].read;
			return true;
		}
		scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapIndex].rightConsensusSequence = GetConsensusSequenceFromFile(consensusFile);

	}
	return true;

}*/



bool RunRacon(char* overlapFile, char* targetFile, char* sequenceFile, char* consensusFile) {
	char command[3000];
	sprintf(command, "racon -w 100 -q 10 -e 0.4 -f %s %s %s > %s 2>log.txt", sequenceFile, overlapFile, targetFile, consensusFile);
	system(command);
	
}

char* GetConsensusSequenceFromFile(char* file) {
	if (access(file, F_OK) != 0) {
		return NULL;
	}

	if (FileIsNull(file) == 1) {
		return NULL;
	}
	FILE* fp;
	if ((fp = fopen(file, "r")) == NULL) {
		printf("%s, does not exist!", file);
		exit(0);
	}

	long int maxSize = 20000;
	char line[20000];
	char* sequence = NULL;
	long int length = 0;
	while ((fgets(line, maxSize, fp)) != NULL) {
		if (line[0] == '>') {
			continue;
		}
		length = strlen(line);
		if (line[length - 1] == '\r' || line[length - 1] == '\n') {
			length--;
		}
		sequence = (char*)malloc(sizeof(char) * (length + 1));
		strncpy(sequence, line, length);
		sequence[length] = '\0';
	}

	fclose(fp);
	return sequence;

}
#endif

