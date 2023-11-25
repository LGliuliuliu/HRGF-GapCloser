#ifndef CONTIGSET_CPP_INCLUDED 
#define CONTIGSET_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "contigSet.h"

using namespace std;

ContigSetHead* GetContigSetFromContigSetFile(char* contigSetFile) {
	ContigSetHead* contigSetHead = (ContigSetHead*)malloc(sizeof(ContigSetHead));
	contigSetHead->contigSet = NULL;
	contigSetHead->contigCount = 0;

	long int maxSize = 90000;
	char* contig = NULL;
	if (NULL == (contig = (char*)malloc(sizeof(char) * maxSize))) {
		perror("contig malloc error");
		exit(1);
	}
	FILE* fp;
	if ((fp = fopen(contigSetFile, "r")) == NULL) {
		printf("%s,does not exist!", contigSetFile);
		exit(0);
	}
	while ((fgets(contig, maxSize, fp)) != NULL) {
		if (contig[0] == '>') {
			contigSetHead->contigCount++;
		}
	}
	fclose(fp);
	contigSetHead->contigSet = (ContigSet*)malloc(sizeof(ContigSet) * contigSetHead->contigCount);
	for (long int i = 0; i < contigSetHead->contigCount; i++) {
		contigSetHead->contigSet[i].contig = NULL;
		contigSetHead->contigSet[i].contigName = NULL;
		contigSetHead->contigSet[i].contigLength = 0;
	}
	if (NULL == (fp = fopen(contigSetFile, "r"))) {
		printf("%s,does not exits!", contigSetFile);
		exit(0);
	}
	long int allocateLength = 0;
	long int contigIndex = -1;
	
	while ((fgets(contig, maxSize, fp)) != NULL) {
		if (contig[0] == '>') {
			if (strlen(contig) == maxSize - 1) {
				while ((fgets(contig, maxSize, fp)) != NULL) {
					if (strlen(contig) != maxSize - 1) {
						break;
					}
				}
			}
			contigIndex++;
			int len = strlen(contig);
			contigSetHead->contigSet[contigIndex].contigName = (char*)malloc(sizeof(char) * len);
			strncpy(contigSetHead->contigSet[contigIndex].contigName, contig + 1, len - 1);
			contigSetHead->contigSet[contigIndex].contigName[len - 2] = '\0';
			continue;
		}

		long int extendLength = strlen(contig);
		if (contig[extendLength - 1] == '\n') {
			extendLength--;
		}
		long int contigLength = 0;
		char* tempContig = NULL;
		if (contigSetHead->contigSet[contigIndex].contig != NULL) {
			if (contigSetHead->contigSet[contigIndex].contigLength + extendLength >= allocateLength) {
				contigLength = contigSetHead->contigSet[contigIndex].contigLength;
				contigSetHead->contigSet[contigIndex].contig = (char*)realloc(contigSetHead->contigSet[contigIndex].contig, allocateLength + maxSize + 1);
				allocateLength = allocateLength + maxSize + 1;
				strncpy(contigSetHead->contigSet[contigIndex].contig + contigSetHead->contigSet[contigIndex].contigLength, contig, extendLength);
				contigSetHead->contigSet[contigIndex].contig[contigLength + extendLength] = '\0';
				contigSetHead->contigSet[contigIndex].contigLength = contigLength + extendLength;

			}
			else {
				strncpy(contigSetHead->contigSet[contigIndex].contig + contigSetHead->contigSet[contigIndex].contigLength, contig, extendLength);
				contigSetHead->contigSet[contigIndex].contig[contigSetHead->contigSet[contigIndex].contigLength + extendLength] = '\0';
				contigSetHead->contigSet[contigIndex].contigLength = contigSetHead->contigSet[contigIndex].contigLength + extendLength;

			}

		}
		else {
			contigSetHead->contigSet[contigIndex].contig = (char*)malloc(sizeof(char) * (maxSize + 1));
			strncpy(contigSetHead->contigSet[contigIndex].contig, contig, extendLength);
			contigSetHead->contigSet[contigIndex].contig[extendLength] = '\0';
			contigSetHead->contigSet[contigIndex].contigLength = extendLength;
			allocateLength = maxSize + 1;

		}


	}

	/*for (int i = 0; i < contigSetHead->contigCount; i++) {
		cout << contigSetHead->contigSet[i].contigName << endl;
		cout << contigSetHead->contigSet[i].contig << endl;
	}*/

	fflush(fp);
	fclose(fp);

	return contigSetHead;
}

//·´Ïò»¥²¹
char* ReverseComplement(char* temp) {
	int len = strlen(temp);
	char* rcTemp = (char*)malloc(sizeof(char) * (len + 1));
	for (int i = 0; i < len; i++) {
		if (temp[i] == 'A' || temp[i] == 'a') {
			rcTemp[len - 1 - i] = 'T';
		}
		else if (temp[i] == 'T' || temp[i] == 't') {
			rcTemp[len - 1 - i] = 'A';
		}
		else if (temp[i] == 'G' || temp[i] == 'g') {
			rcTemp[len - 1 - i] = 'C';
		}
		else if (temp[i] == 'C' || temp[i] == 'c') {
			rcTemp[len - 1 - i] = 'G';
		}
		else if (temp[i] == 'N' || temp[i] == 'n') {
			rcTemp[len - 1 - i] = 'N';
		}
	}
	rcTemp[len] = '\0';
	return rcTemp;
}

void ReverseSequence(char* temp) {
	int len = strlen(temp);
	char a;
	for (int i = 0; i < len / 2; i++) {
		a = temp[i];
		temp[i] = temp[len - i - 1];
		temp[len - i - 1] = a;
	}

}


long int GetMinValue(long int a, long int b) {
	if (a > b) {
		return b;
	}
	return a;
}













#endif