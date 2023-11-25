
#ifndef SCAFFOLDSET_CPP_INCLUDED 
#define SCAFFOLDSET_CPP_INCLUDED 

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


#include"scaffoldSet.h"

using namespace std;

void DelectShortContig(char* scaffoldSetFile, char* longScaffoldSetFile) {
	ScaffoldSetHead* scaffoldSetHead = GetScaffoldSetFromScaffoldFile(scaffoldSetFile);//获取scaffold以及GAP位置信息
	GetContigSetFromScaffoldSetHead(scaffoldSetHead);
	FILE* fp;
	if ((fp = fopen(longScaffoldSetFile, "w")) == NULL) {
		printf("%s, does not exist!", longScaffoldSetFile);
		exit(0);
	}
	int flag = 0;
	long int length = 1000000;
	char* contig = (char*)malloc(sizeof(char) * length);
	long int cutReadCount = 200;
	long int minReadLength = 6000;

	/*for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {
			flag = 1;
		}
		cout << "i=" << i << endl;
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {
			cout << "j=" << j << endl;
			cout <<"contig:" << scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate << "--------" << scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigEndCoordinate << "-------" << scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength << endl;
			if (j + flag < scaffoldSetHead->scaffoldSet[i].gapCount) {
				cout <<"gap:" << scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapStartCoordinate << "-------" << scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapEndCoordinate << "--------" << scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength << endl;
			}
		}
		cout << endl;
	}*/

	//删除两端读数
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {
			flag = 1;
		}

		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {

			if (j == 0 || j == scaffoldSetHead->scaffoldSet[i].contigCount - 1) {
				if (j == 0 && scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength > minReadLength) {
					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigEndCoordinate = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigEndCoordinate - cutReadCount;
					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength - cutReadCount;

					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapStartCoordinate = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapStartCoordinate - cutReadCount;
					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength + cutReadCount;
				}
				if (j == scaffoldSetHead->scaffoldSet[i].contigCount - 1 && scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength > minReadLength) {
					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate + cutReadCount;
					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength - cutReadCount;

					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapEndCoordinate = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapEndCoordinate + cutReadCount;
					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapLength = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapLength + cutReadCount;

				}

			}
			else {
				if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength > minReadLength) {
					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate + cutReadCount;
					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength - cutReadCount;

					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapEndCoordinate = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapEndCoordinate + cutReadCount;
					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapLength = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag - 1].gapLength + cutReadCount;


					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigEndCoordinate = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigEndCoordinate - cutReadCount;
					scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength = scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength - cutReadCount;

					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapStartCoordinate = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapStartCoordinate - cutReadCount;
					scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength = scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength + cutReadCount;
				}
			}
		}
	}

	flag = 0;
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			fprintf(fp, ">%ld\n%s\n", i, scaffoldSetHead->scaffoldSet[i].scaffold);
			continue;
		}
		fprintf(fp, ">%ld\n", i);

		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {
			long int t = 0;
			if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapLength < length - 1) {
				for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[0].gapLength; t++) {
					contig[t] = 'N';
				}
				contig[t] = '\0';
				fprintf(fp, "%s", contig);
			}
			else {
				char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].gapSet[0].gapLength + 1);
				for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[0].gapLength; t++) {
					tempContig[t] = 'N';
				}
				tempContig[t] = '\0';
				fprintf(fp, "%s", tempContig);
			}
			flag = 1;
			
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {
			
			if (j == 0 || j== scaffoldSetHead->scaffoldSet[i].contigCount-1) {//输出contig
				if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength < length - 1) {
					strncpy(contig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
					contig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
					fprintf(fp, "%s", contig);
				}
				else {
					char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength + 1);
					strncpy(tempContig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
					tempContig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
					fprintf(fp, "%s", tempContig);
				}
				
			}
			else {
				if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength > 1000) {//删除长度短与2000的contig
					if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength < length - 1) {
						strncpy(contig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
						contig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
						fprintf(fp, "%s", contig);
					}
					else {
						char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength + 1);
						strncpy(tempContig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
						tempContig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
						fprintf(fp, "%s", tempContig);
					}
				}
				else {
					long int t = 0;
					for (t = 0; t < scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength; t++) {
						contig[t] = 'N';
					}
					contig[t] = '\0';
					fprintf(fp, "%s", contig);
				}
				
			}
			if (j + flag < scaffoldSetHead->scaffoldSet[i].gapCount) {
				long int t = 0;
				if (scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength < length - 1) {
					for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength; t++) {
						contig[t] = 'N';
					}
					contig[t] = '\0';
					fprintf(fp, "%s", contig);
				}
				else {
					char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength + 1);
					for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength; t++) {
						tempContig[t] = 'N';
					}
					tempContig[t] = '\0';
					fprintf(fp, "%s", tempContig);
				}


			}
		}
		fprintf(fp, "\n");
	}
	fflush(fp);
	fclose(fp);

}

//获取每个scaffold的gap信息及数量
GapSet* GetGapSetFromOneScaffold(char* scaffold, long int& gapCount) {
	long int scaffoldLength = strlen(scaffold);
	gapCount = 0;
	bool token = true;
	for (long int i = 0; i < scaffoldLength; i++) {
		if ((scaffold[i] == 'N' || scaffold[i] == 'n') && token != false) {
			gapCount++;
			token = false;
		}
		if (scaffold[i] != 'N' && scaffold[i] != 'n' && token != true) {
			token = true;
		}
	}
	GapSet* gapSet = (GapSet*)malloc(sizeof(GapSet) * gapCount);
	for (long int i = 0; i < gapCount; i++) {
		gapSet[i].gapStartCoordinate = -1;
		gapSet[i].gapEndCoordinate = -1;
		gapSet[i].gapLength = -1;
		gapSet[i].estimatedGapDistance = 0;
		gapSet[i].estimatedReadCount = 0;
	}
	gapCount = 0;
	token = true;
	for (long int i = 0; i < scaffoldLength; i++) {
		if ((scaffold[i] == 'N' || scaffold[i] == 'n') && token != false) {
			gapSet[gapCount].gapStartCoordinate = i;
			gapCount++;
			token = false;
		}
		if (scaffold[i] != 'N' && scaffold[i] != 'n' && token != true) {
			if (gapCount > 0) {
				gapSet[gapCount - 1].gapEndCoordinate = i - 1;
			}
			token = true;
		}
	}
	if (token == false) {
		gapSet[gapCount - 1].gapEndCoordinate = scaffoldLength - 1;
	}
	for (long int i = 0; i < gapCount; i++) {
		gapSet[i].gapLength = gapSet[i].gapEndCoordinate - gapSet[i].gapStartCoordinate + 1;
	}
	/*for (int i = 0; i < gapCount; i++) {
		cout << i << "-" << gapCount << endl;
		cout << gapSet[i].gapStartCoordinate << "+" << gapSet[i].gapEndCoordinate << "---" << gapSet[i].gapLength<< endl;
	}*/
	return gapSet;
}

//获取scaffold名称、序列
ScaffoldSetHead* GetScaffoldSetFromScaffoldFile(char* scaffoldSetFile) {
	long int i = 0;
	long int j = 0;
	long int maxSize = 1000000;
	char* scaffold = NULL;
	if (NULL == (scaffold = (char*)malloc(sizeof(char)*maxSize))) {
		perror("scaffold malloc error");
		exit(1);
	}
	long int scaffoldCount = 0;
	FILE* fp;
	if ((fp = fopen(scaffoldSetFile, "r")) == NULL) {
		printf("%s,does not exist!", scaffoldSetFile);
		exit(0);
	}
	while ((fgets(scaffold, maxSize, fp)) != NULL) {
		if (scaffold[0] == '>') {
			scaffoldCount++;
		}
	}
	
	fclose(fp);
	ScaffoldSetHead* scaffoldSetHead = NULL;
	if (NULL == (scaffoldSetHead = (ScaffoldSetHead*)malloc(sizeof(ScaffoldSetHead)))) {
		perror("scaffoldSetHead malloc error!");
		exit(1);
	}
	scaffoldSetHead->scaffoldCount = scaffoldCount;
	scaffoldSetHead->allGapCount = 0;
	if (NULL == (scaffoldSetHead->scaffoldSet = (ScaffoldSet*)malloc(sizeof(ScaffoldSet)*scaffoldCount))) {
		perror("scaffoldSet malloc error!");
		exit(1);
	}
	
	for (i = 0; i < scaffoldCount; i++) {
		scaffoldSetHead->scaffoldSet[i].scaffold = NULL;
		scaffoldSetHead->scaffoldSet[i].scaffoldName = NULL;
		scaffoldSetHead->scaffoldSet[i].scaffoldLength = 0;
		scaffoldSetHead->scaffoldSet[i].gapSet = NULL;
		scaffoldSetHead->scaffoldSet[i].contigPostion = NULL;
		scaffoldSetHead->scaffoldSet[i].gapCount = 0;
		scaffoldSetHead->scaffoldSet[i].contigCount = 0;
		scaffoldSetHead->scaffoldSet[i].hifiGapSet = NULL;
	}
	if ((fp = fopen(scaffoldSetFile, "r")) == NULL) {
		printf("%s,does not exist!", scaffoldSetFile);
		exit(0);
	}
	long int scaffoldIndex = -1;
	long int allocateLength = 0;
	while ((fgets(scaffold, maxSize, fp)) != NULL) {//最多读取maxSize-1个字符
		if (scaffold[0] == '>') {
			if (strlen(scaffold) == maxSize - 1) {
				while ((fgets(scaffold, maxSize, fp)) != NULL) {
					if (strlen(scaffold) != maxSize - 1) {
						break;
					}
					
				}
			}
			if (scaffoldIndex >= 0) {
				scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet = GetGapSetFromOneScaffold(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, scaffoldSetHead->scaffoldSet[scaffoldIndex].gapCount);//获取gap对应信息
			}
			scaffoldIndex++;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName = (char*)malloc(sizeof(char) * 20);
			sprintf(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName, ">scaffold_%ld", scaffoldIndex);
			continue;
		}
		long int extendLength = strlen(scaffold);
		if (scaffold[extendLength - 1] == '\n') {
			extendLength--;
		}
		long int scaffoldLength = 0;
		char* tempScaffold = NULL;
		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold != NULL) {
			if (scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength + extendLength >= allocateLength) {
				scaffoldLength = scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength;
				tempScaffold = (char*)malloc(sizeof(char) * (scaffoldLength + 1));
				strncpy(tempScaffold, scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, scaffoldLength);
				free(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold);
				scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold = (char*)malloc(sizeof(char) * (allocateLength + maxSize + 1));
				strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, tempScaffold, scaffoldLength);
				allocateLength = allocateLength + maxSize + 1;
				strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold + scaffoldLength, scaffold, extendLength);
				scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[scaffoldLength + extendLength] = '\0';
				scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = scaffoldLength + extendLength;
				free(tempScaffold);

			}
			else {
				strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold + scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength, scaffold, extendLength);
				scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength + extendLength] = '\0';
				scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength + extendLength;

			}

		}
		else {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold = (char*)malloc(sizeof(char) * (maxSize + 1));
			strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, scaffold, extendLength);
			scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[extendLength] = '\0';
			scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = extendLength;
			allocateLength = maxSize + 1;
		}


	}
	scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet = GetGapSetFromOneScaffold(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, scaffoldSetHead->scaffoldSet[scaffoldIndex].gapCount);
	fclose(fp);
	//for (int k = 0; k < scaffoldCount; k++) {//输出gap在scaffold上开始与结束的位置信息
	//	
	//	for (int j = 0; j < scaffoldSetHead->scaffoldSet[k].gapCount; j++) {
	//		cout << "k=" << k <<"     j="<<j<<"   gapStart="<< scaffoldSetHead->scaffoldSet[k].gapSet[j].gapStartCoordinate<<"  gapend=" << scaffoldSetHead->scaffoldSet[k].gapSet[j].gapEndCoordinate<< "       gapLength=" << scaffoldSetHead->scaffoldSet[k].gapSet[j].gapLength << endl;
	//		//cout << scaffoldSetHead->scaffoldSet[k].gapSet[j].gapStartCoordinate<<"------"<< scaffoldSetHead->scaffoldSet[k].gapSet[j].gapEndCoordinate<< endl;
	//	}
	//	cout << endl;
	//}
	//cout << scaffoldSetHead->scaffoldSet[6664].scaffoldName << endl;
	//cout << scaffoldSetHead->scaffoldSet[6664].scaffold << endl;
    return scaffoldSetHead;
}

//获取contig在scaffold中的位置信息
void GetContigSetFromScaffoldSetHead(ScaffoldSetHead* scaffoldSetHead) {
	long int gapCount = 0;
	long int contigCount = 0;
	long int scaffoldLength = 0;
	GapSet* gapSet = NULL;
	for (int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		gapCount = scaffoldSetHead->scaffoldSet[i].gapCount;
		gapSet = scaffoldSetHead->scaffoldSet[i].gapSet;
		scaffoldLength = scaffoldSetHead->scaffoldSet[i].scaffoldLength;
		contigCount = gapCount + 1;
		long int start = 0;
		long int end = 0;
		if (gapCount > 0 && gapSet[0].gapStartCoordinate == 0) {
			contigCount--;
			start = 1;
		}
		if (gapCount > 0 && gapSet[gapCount - 1].gapEndCoordinate == scaffoldLength - 1) {
			contigCount--;
			end = 1;
		}
		if (contigCount < 1) {//没有gap的情况
			continue;
		}
		ContigPosition* contigPosition = (ContigPosition*)malloc(sizeof(ContigPosition)*contigCount);
		for (long int j = 0; j < contigCount; j++) {
			contigPosition[j].contigStartCoordinate = -1;
			contigPosition[j].contigEndCoordinate = -1;
			contigPosition[j].contigLength = -1;
		}
		for (long int j = 0; j < gapCount; j++) {
			if (gapSet[j].gapStartCoordinate != 0 && gapSet[j].gapEndCoordinate != scaffoldLength - 1) {
				contigPosition[j - start].contigEndCoordinate = gapSet[j].gapStartCoordinate - 1;
				contigPosition[j - start + 1].contigStartCoordinate = gapSet[j].gapEndCoordinate + 1;
			}
			else if (gapSet[j].gapStartCoordinate == 0) {
				contigPosition[j].contigStartCoordinate = gapSet[j].gapEndCoordinate + 1;
			}
			else {
				contigPosition[j - start].contigEndCoordinate = gapSet[j].gapStartCoordinate - 1;
			}
		}
		if (start == 0) {
			contigPosition[0].contigStartCoordinate = 0;
		}
		if (end == 0) {
			contigPosition[contigCount - 1].contigEndCoordinate = scaffoldLength - 1;
		}
		for (long int j = 0; j < contigCount; j++) {
			contigPosition[j].contigLength = contigPosition[j].contigEndCoordinate - contigPosition[j].contigStartCoordinate + 1;
		}
		scaffoldSetHead->scaffoldSet[i].contigPostion = contigPosition;
		scaffoldSetHead->scaffoldSet[i].contigCount = contigCount;

	}
	//long int num = 0;
	//for (int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
	//	if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
	//		continue;
	//	}
	//	for (int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {
	//		cout << "i=" << i << "    j=" << j <<"  num=" <<num<< "    contigLenth=" << scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength << endl;
	//			/*if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength <= 4000) {
	//			cout <<"j="<<j<<"    contig:"<< scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate << "-----" << scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigEndCoordinate << "=====" << scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength << endl;
	//			cout << endl;
	//		}*/
	//		num++;
	//	}
	//	cout << endl;
	//}
	
}

void OutPutContigSetInScaffoldSetHead(ScaffoldSetHead* scaffoldSetHead, char* contigSetFile) {
	FILE* fp;
	if ((fp = fopen(contigSetFile, "w")) == NULL) {
		printf("%s,does exist!", contigSetFile);
		exit(0);
	}
	long int contigCount = 0;
	long int length = 3000000;
	char* contig = (char*)malloc(sizeof(char) * length);
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {//注意       没有gap的不保存
			continue;
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {
			fprintf(fp, ">contig%ld_%ld-%ld\n",contigCount, i, j );
			contigCount++;
			if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength < length - 1) {
				strncpy(contig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				contig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s\n", contig);
			}else {
				char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength + 1);
				strncpy(tempContig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				tempContig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s\n", tempContig);
			}
		}
	}
	fflush(fp);
	fclose(fp);


	//for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
	//	if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {//注意       没有gap的不保存
	//		continue;
	//	}
	//	cout << "i=" << i << "----" << scaffoldSetHead->scaffoldSet[i].contigCount << endl;
	//}
	

}

//通过比对信息获取gap的平均长度以及支持该gap的读数的个数
void GetEstimatedGapLength(ScaffoldSetHead* scaffoldSetHead, ContigSetHead* contigSetHead, char* file, char* line, int maxSize) {//file:(contigPathInHifiRead.fa)
	FILE* fp;
	if ((fp = fopen(file, "r")) == NULL) {
		printf("%s,does not exist!", file);
		exit(0);
	}
	const char* split = ",";
	char* p;
	long int dataNum = 7;//数据数量
	int maxCount = 100;
	int aligningSingle[dataNum * maxCount];
	long int readIndex=-1;
	long int scaffoldCount = 6640;
	long int gapCount = 11;
	while ((fgets(line, maxSize, fp)) != NULL) {
	/*for(int k=0;k<2;k++){
		fgets(line, maxSize, fp);*/
		p = strtok(line, split);
		int count = atoi(p);
		if (count <= 1) {//********************************************改
			continue;
		}
		p = strtok(NULL, split);
		readIndex = atoi(p);
		p = strtok(NULL, split);
		int readLength = atoi(p);
		int a = 1;
		int b = 0;
		
		while (a <= count * dataNum) {
			p = strtok(NULL, split);
			aligningSingle[b] = atoi(p);
			b++;
			a++;
		}
		/*cout << "count=" << count << endl;
		for (long int i = 0; i < count * dataNum; i++) {
			cout << aligningSingle[i] << ",";
		}
		cout << endl;
		cout << endl;*/
		
		long int scaffoldIndex = -1;
		long int gapInScaffoldIndex = -1;
		long int gapIndex = -1;
		if (count == 1) {
			continue;
		}
		int gapDistance = 0;

		for (int i = 0; i < count - 1; i++) {
			for (int j = i + 1; j < count; j++) {
				scaffoldIndex = -1;
				gapIndex = -1;
				if (abs(aligningSingle[i * dataNum] - aligningSingle[j * dataNum]) == 1 && aligningSingle[i * dataNum + 6] == aligningSingle[j * dataNum + 6] 
					&& aligningSingle[i * dataNum + 6] == 1 && aligningSingle[i * dataNum + 1] < aligningSingle[j * dataNum + 1]) {
					long int contigLength = contigSetHead->contigSet[aligningSingle[i * dataNum]].contigLength;
					if (contigLength - aligningSingle[i * dataNum + 4] - 1 < 8000 && aligningSingle[j * dataNum + 3] < 8000) {
						gapIndex = GetGapIndexFromContigIndex(scaffoldSetHead, aligningSingle[i * dataNum], aligningSingle[j * dataNum], scaffoldIndex, gapInScaffoldIndex);
					}
				}
				if (abs(aligningSingle[i * dataNum] - aligningSingle[j * dataNum]) == 1 && aligningSingle[i * dataNum + 6] == aligningSingle[j * dataNum + 6] 
					&& aligningSingle[i * dataNum + 6] == 0  && aligningSingle[i * dataNum + 1] > aligningSingle[j * dataNum + 1]) {
					long int contigLength = contigSetHead->contigSet[aligningSingle[i * dataNum]].contigLength;
					if (contigLength - aligningSingle[i * dataNum + 4] - 1 < 8000 && aligningSingle[j * dataNum + 3] < 8000) {
						gapIndex = GetGapIndexFromContigIndex(scaffoldSetHead, aligningSingle[i * dataNum], aligningSingle[j * dataNum], scaffoldIndex, gapInScaffoldIndex);
					}
				}
				long int diff = 200000;
				if (gapIndex != -1) {
					if (aligningSingle[i * dataNum + 6] == 1) {
						diff = abs(scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].gapLength - abs(aligningSingle[j * dataNum + 1] - aligningSingle[i * dataNum + 2]) - 1);
					}
					else {
						diff = abs(scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].gapLength - abs(aligningSingle[i * dataNum + 1] - aligningSingle[j * dataNum + 2]) - 1);
					}
				}
				if (gapIndex != -1 && diff < 10000) {//************************************************************同向且gap误差小于700
					if (aligningSingle[i * dataNum + 6] == 1) {
						scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance
						+ abs(aligningSingle[j * dataNum + 1] - aligningSingle[i * dataNum + 2]) + 1;

					}
					else {
						scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance
							+ abs(aligningSingle[i * dataNum + 1] - aligningSingle[j * dataNum + 2]) + 1;
					}
					
					scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedReadCount++;//估计的读数数量
				}
				
			}
		}
	}

	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {//求gap平均值
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
			if (scaffoldSetHead->scaffoldSet[i].gapSet[j].estimatedReadCount > 1) {
				scaffoldSetHead->scaffoldSet[i].gapSet[j].estimatedGapDistance = scaffoldSetHead->scaffoldSet[i].gapSet[j].estimatedGapDistance / scaffoldSetHead->scaffoldSet[i].gapSet[j].estimatedReadCount;//平均gap长度
			}else {
				scaffoldSetHead->scaffoldSet[i].gapSet[j].estimatedGapDistance = scaffoldSetHead->scaffoldSet[i].gapSet[j].gapLength;
			}
		}
	}
	fclose(fp);
}

//获得比对到gap区域读数的位置
void GetReadPoinstionInGap(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, char* file, char* line, int maxSize, FILE* hifiFileFP) {//file:contigPathInHifiRead.fa

	long int extendLength = 0;
	long int diff = 200;
	long int extendReadLength = 10;
	
	FILE* fp;
	if ((fp = fopen(file, "r")) == NULL) {
		printf("%s, does not exist!", file);
		exit(0);
	}
	const char* split = ",";
	char* p;

	long int dataNum = 7;
	int maxCount = 100;
	long int aligningSingle[dataNum * maxCount];
	long int aligningSingle01[dataNum * maxCount];
	long int readIndex = -1;
	while ((fgets(line, maxSize, fp)) != NULL) {
	/*for(int k=0;k<50;k++){
		fgets(line, maxSize, fp);*/
		p = strtok(line, split);
		int count = atoi(p);
		if (count < 1) {
			continue;
		}
		p = strtok(NULL, split);
		readIndex = atoi(p);

		p = strtok(NULL, split);
		int readLength = atoi(p);

		//cout << "count=" << count << "  readIndex=" << readIndex << " readLength=" << readLength << endl;
		int a = 1;
		int b = 0;
		while (a <= count * dataNum) {
			p = strtok(NULL, split);
			aligningSingle[b] = atoi(p);
			b++;
			a++;
		}
		/*for (int i = 0; i < count * dataNum; i++) {
			cout << aligningSingle[i] << ",";
		}
		cout << endl;
		cout << endl;*/


		long int j = 0;
		long int num = 0;

		long int startFlag = 0;
		long int endFlag = 0;

		long int scaffoldIndex = -1;
		long int gapInScaffoldIndex = -1;
		long int gapIndex = -1;

		long int leftContigDiff = -1;
		long int rightContigDiff = -1;
		long int leftContigOverlap = -1;
		long int rightContigOverlap = -1;
		long int weigth = 0;
		//单一contig获取左右读数
		if (count == 1) {//只有一个读数
			
			gapIndex = -1;
			gapInScaffoldIndex = -1;
			scaffoldIndex = -1;
			leftContigDiff = -1;
			rightContigDiff = -1;
			leftContigOverlap = -1;
			rightContigOverlap = -1;
			weigth = 1;
			
			gapIndex = GetGapIndexFromContigIndexOne(scaffoldSetHead, aligningSingle[0], 0, scaffoldIndex, gapInScaffoldIndex);//计算contig左边gap索引-正向比对
			//cout << "RscaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
			if (gapIndex != -1) {
				int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
				long int contigLength = contigSetHead->contigSet[aligningSingle[0]].contigLength;
				//cout << "RreadLength=" << readLength << "       contigLength=" << contigLength << endl;
				//cout << "RstartDiff=" << aligningSingle[3]<< endl;
				if (aligningSingle[6] == 1 && aligningSingle[3] < aligningSingle[1] && aligningSingle[3] < 600 && aligningSingle[1]>100) {

					//cout << "R1*************************************************************contigIndex=" << aligningSingle[0] << endl;
					int start = 0;
					int end = 0;
					rightContigDiff = aligningSingle[3];
					rightContigOverlap = aligningSingle[5];
					if (aligningSingle[1] - 1 > distance + diff) {
						start = aligningSingle[1] - distance - diff;
					}
					else {
						start = 0;
					}
					end = aligningSingle[1] - 1;
					fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle[6]);
				}

				if (aligningSingle[6] == 0 && aligningSingle[3] < readLength - aligningSingle[2] - 1 && aligningSingle[3] < 600 && readLength - aligningSingle[2] - 1>100) {//计算contig左边gap索引-反向比对

					//cout << "R2************************************************************************contigIndex=" << aligningSingle[0] << endl;
					int start = 0;
					int end = 0;
					rightContigDiff = aligningSingle[3];
					rightContigOverlap = aligningSingle[5];
					if (readLength - aligningSingle[2] - 1 > distance + diff) {
						end = aligningSingle[2] + distance + diff;
					}
					else {
						end = readLength - 1;
					}
					start = aligningSingle[2] + 1;
					fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle[6]);

				}
				//cout << endl;

			}
			
			gapIndex = -1;
			gapInScaffoldIndex = -1;
			scaffoldIndex = -1;
			leftContigDiff = -1;
			rightContigDiff = -1;
			leftContigOverlap = -1;
			rightContigOverlap = -1;
			gapIndex = -1;
			gapIndex = GetGapIndexFromContigIndexOne(scaffoldSetHead, aligningSingle[0], 1, scaffoldIndex, gapInScaffoldIndex);//计算contig右边gap索引-正向比对
			//cout << "LscaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
			
			if (gapIndex != -1) {
				int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
				long int contigLength = contigSetHead->contigSet[aligningSingle[0]].contigLength;
				//cout << "LreadLength=" << readLength << "       contigLength=" << contigLength << endl;
				//cout << "LstartDiff=" << contigLength - aligningSingle[4] - 1 << endl;
				
				if (aligningSingle[6] == 1 && contigLength - aligningSingle[4] - 1 < readLength - aligningSingle[2] - 1 && contigLength - aligningSingle[4] - 1 < 600 && readLength - aligningSingle[2] - 1>100) {
					//cout << "L1************************************************************************contigIndex=" << aligningSingle[0] << endl;
					int start = 0;
					int end = 0;
					leftContigDiff = contigLength - aligningSingle[4] - 1;
					leftContigOverlap = aligningSingle[5];
					if (readLength - aligningSingle[2] - 1 > distance + diff) {
						end = aligningSingle[2] + distance + diff;

					}
					else {
						end = readLength - 1;
					}
					start = aligningSingle[2] + 1;
					fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle[6]);
				}
				if (aligningSingle[6] == 0 && contigLength - aligningSingle[4] - 1 < aligningSingle[1] && contigLength - aligningSingle[4] - 1 < 600 && aligningSingle[1]>100) {
					//cout << "L2************************************************************************contigIndex=" << aligningSingle[0] << endl;
					int start = 0;
					int end = 0;
					leftContigDiff = contigLength - aligningSingle[4] - 1;
					leftContigOverlap = aligningSingle[5];
					if (aligningSingle[1] > distance + diff) {
						start = aligningSingle[1] - distance - diff;
					}
					else {
						start = 0;
					}
					end = aligningSingle[1] - 1;
					fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle[6]);
				}
				//cout << endl;
			}
			continue;
		}

		for (int i = 0; i < count; i++) {

			if (abs(aligningSingle[i * dataNum] - aligningSingle[(i + 1) * dataNum])==1 && i!=count-1) {
				aligningSingle01[num * dataNum] = aligningSingle[i * dataNum];
				aligningSingle01[num * dataNum+1] = aligningSingle[i * dataNum+1];
				aligningSingle01[num * dataNum+2] = aligningSingle[i * dataNum+2];
				aligningSingle01[num * dataNum+3] = aligningSingle[i * dataNum+3];
				aligningSingle01[num * dataNum+4] = aligningSingle[i * dataNum+4];
				aligningSingle01[num * dataNum+5] = aligningSingle[i * dataNum+5];
				aligningSingle01[num * dataNum+6] = aligningSingle[i * dataNum+6];
				num++;
				continue;
			}
			if (abs(aligningSingle[i * dataNum] - aligningSingle[(i + 1) * dataNum]) == 0 && i != count - 1) {
				aligningSingle01[num * dataNum] = aligningSingle[i * dataNum];
				aligningSingle01[num * dataNum + 1] = aligningSingle[i * dataNum + 1];
				aligningSingle01[num * dataNum + 2] = aligningSingle[i * dataNum + 2];
				aligningSingle01[num * dataNum + 3] = aligningSingle[i * dataNum + 3];
				aligningSingle01[num * dataNum + 4] = aligningSingle[i * dataNum + 4];
				aligningSingle01[num * dataNum + 5] = aligningSingle[i * dataNum + 5];
				aligningSingle01[num * dataNum + 6] = aligningSingle[i * dataNum + 6];
				num++;
				continue;
			}
			if (abs(aligningSingle[i * dataNum] - aligningSingle[(i + 1) * dataNum]) > 1 || i==count-1) {
				aligningSingle01[num * dataNum] = aligningSingle[i * dataNum];
				aligningSingle01[num * dataNum + 1] = aligningSingle[i * dataNum + 1];
				aligningSingle01[num * dataNum + 2] = aligningSingle[i * dataNum + 2];
				aligningSingle01[num * dataNum + 3] = aligningSingle[i * dataNum + 3];
				aligningSingle01[num * dataNum + 4] = aligningSingle[i * dataNum + 4];
				aligningSingle01[num * dataNum + 5] = aligningSingle[i * dataNum + 5];
				aligningSingle01[num * dataNum + 6] = aligningSingle[i * dataNum + 6];
				num++;
			}

			/*for (int a = 0; a < num * dataNum; a++) {
				cout << aligningSingle01[a] << ",";
			}

			cout << "num=" << num << endl;*/
			weigth = num;
			long int scaffoldIndexNum = 6635;
			//许多contig中的单独contig获取左右contig
			if (num ==1 ) {
				
				gapIndex = -1;
				gapInScaffoldIndex = -1;
				scaffoldIndex = -1;
				leftContigDiff = -1;
				rightContigDiff = -1;
				leftContigOverlap = -1;
				rightContigOverlap = -1;
				weigth = num;
				
				gapIndex = GetGapIndexFromContigIndexOne(scaffoldSetHead, aligningSingle01[0], 0, scaffoldIndex, gapInScaffoldIndex);//计算contig左边gap索引-正向比对
				//cout << "LscaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
				if (gapIndex != -1) {
					
					int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
					long int contigLength = contigSetHead->contigSet[aligningSingle01[0]].contigLength;
					//cout << "LreadLength=" << readLength << "       contigLength=" << contigLength << endl;
					//cout << "LstartDiff=" << aligningSingle01[3]<< endl;
					if (aligningSingle01[6] == 1 && aligningSingle01[3]<aligningSingle01[1] && aligningSingle01[3]<600 && aligningSingle01[1]>100) {
						
						//cout <<"L1*************************************************************contigIndex="<< aligningSingle01[0]<< endl;
						int start = 0;
						int end = 0;
						rightContigDiff = aligningSingle01[3];
						rightContigOverlap = aligningSingle01[5];
						if (aligningSingle01[1] - 1 > distance + diff) {
							start = aligningSingle01[1] - distance - diff;
						}
						else {
							start = 0;
						}
						end = aligningSingle01[1] - 1;
						/*if (scaffoldIndex == scaffoldIndexNum) {
							cout << "**********************************************************scaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,r\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
						}*/
						fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
					}

					if (aligningSingle01[6] == 0 && aligningSingle01[3] < readLength-aligningSingle01[2]-1 && aligningSingle01[3] < 600 && readLength - aligningSingle01[2] - 1>100) {//计算contig左边gap索引-反向比对
						
						//cout << "L2************************************************************************contigIndex=" << aligningSingle01[0] << endl;
						int start = 0;
						int end = 0;
						rightContigDiff = aligningSingle01[3];
						rightContigOverlap = aligningSingle01[5];
						if (readLength - aligningSingle01[2] - 1 > distance + diff) {
							end = aligningSingle01[2] + distance + diff;
						}
						else {
							end = readLength-1;
						}
						start = aligningSingle01[2] + 1;
						/*if (scaffoldIndex == scaffoldIndexNum) {
							cout << "**********************************************************scaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,r\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
						}*/
						fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);

					}
					//cout << endl;

				}
				gapIndex = -1;
				gapInScaffoldIndex = -1;
				scaffoldIndex = -1;
				leftContigDiff = -1;
				rightContigDiff = -1;
				leftContigOverlap = -1;
				rightContigOverlap = -1;
				gapIndex = -1;

				gapIndex = GetGapIndexFromContigIndexOne(scaffoldSetHead, aligningSingle01[0], 1, scaffoldIndex, gapInScaffoldIndex);//计算contig右边gap索引-正向比对
				//cout << "LscaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
				if (gapIndex != -1) {
					
					int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
					long int contigLength = contigSetHead->contigSet[aligningSingle01[0]].contigLength;
					//cout << "LreadLength=" << readLength << "       contigLength=" << contigLength << endl;
					//cout << "LstartDiff=" << contigLength - aligningSingle01[4] - 1 << endl;
					if (aligningSingle01[6] == 1 && contigLength - aligningSingle01[4] - 1 < readLength - aligningSingle01[2] - 1 && contigLength - aligningSingle01[4] - 1 < 600 && readLength - aligningSingle01[2] - 1>100) {
						//cout << "L1************************************************************************contigIndex=" << aligningSingle01[0] << endl;
						int start = 0;
						int end = 0;
						leftContigDiff = contigLength - aligningSingle01[4] - 1;
						leftContigOverlap = aligningSingle01[5];
						if (readLength - aligningSingle01[2] - 1 > distance + diff) {
							end = aligningSingle01[2] + distance + diff;

						}
						else {
							end = readLength - 1;
						}
						start = aligningSingle01[2] + 1;
						/*if (scaffoldIndex == scaffoldIndexNum) {
							cout << "**********************************************************scaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,l\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
						}*/
						fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
					}
					if (aligningSingle01[6] == 0 && contigLength - aligningSingle01[4] - 1 < aligningSingle01[1] && contigLength - aligningSingle01[4] - 1 < 600 && aligningSingle01[1]>100) {
						//cout << "L2************************************************************************contigIndex=" << aligningSingle01[0] << endl;
						int start = 0;
						int end = 0;
						leftContigDiff = contigLength - aligningSingle01[4] - 1;
						leftContigOverlap = aligningSingle01[5];
						if (aligningSingle01[1] > distance + diff) {
							start = aligningSingle01[1] - distance - diff;
						}
						else {
							start = 0;
						}
						end = aligningSingle01[1] - 1;
						/*if (scaffoldIndex == scaffoldIndexNum) {
							cout << "**********************************************************scaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,l\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
						}*/
						fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
					}
					//cout << endl;
					
				}
				
				num = 0;
				continue;
			}
			int n = 0;
			for (int m = 0; m < num - 1; m++) {
				n = m + 1;
				//cout << "m=" << m<< "  n=" << n << endl;
				gapIndex = -1;
				gapInScaffoldIndex = -1;
				scaffoldIndex = -1;
				leftContigDiff = -1;
				rightContigDiff = -1;
				leftContigOverlap = -1;
				rightContigOverlap = -1;

				//cout << "contig0=" << aligningSingle01[m * dataNum] << "--------" << "contig0=" << aligningSingle01[n * dataNum] << endl;

				//连续contig获取第一个contig的开始读数，及最后一个contig的最后读数
				if (m == 0) {//计算第一个contig左边gap的对应序列
					
					gapIndex = GetGapIndexFromContigIndexOne(scaffoldSetHead, aligningSingle01[0], 0, scaffoldIndex, gapInScaffoldIndex);
					//cout << "RscaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
					if (gapIndex != -1) {
						int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
						long int contigLength = contigSetHead->contigSet[aligningSingle01[0]].contigLength;
						//cout << "RreadLength=" << readLength << "       contigLength=" << contigLength << "   Diff=" << aligningSingle01[3] << endl;
						if (aligningSingle01[6]== aligningSingle01[13] && aligningSingle01[6] == 1 && aligningSingle01[3] < aligningSingle01[1] && aligningSingle01[3] < 600 && aligningSingle01[1]>100) {

							//cout <<"R1*************************************************************contigIndex="<< aligningSingle01[0]<< endl;
							int start = 0;
							int end = 0;
							rightContigDiff = aligningSingle01[3];
							rightContigOverlap = aligningSingle01[5];
							if (aligningSingle01[1] - 1 > distance + diff) {
								start = aligningSingle01[1] - distance - diff;
							}
							else {
								start = 0;
							}
							end = aligningSingle01[1] - 1;
							/*if (scaffoldIndex == scaffoldIndexNum) {
								fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R0\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
							}*/
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
						}

						if (aligningSingle01[6] == aligningSingle01[13] && aligningSingle01[6] == 0 && aligningSingle01[3] < readLength - aligningSingle01[2] - 1 && aligningSingle01[3] < 600 && readLength - aligningSingle01[2] - 1>100) {//计算contig左边gap索引-反向比对

							//cout << "R2************************************************************************contigIndex=" << aligningSingle01[0] << endl;
							int start = 0;
							int end = 0;
							rightContigDiff = aligningSingle01[3];
							rightContigOverlap = aligningSingle01[5];
							if (readLength - aligningSingle01[2] - 1 > distance + diff) {
								end = aligningSingle01[2] + distance + diff;
							}
							else {
								end = readLength - 1;
							}
							start = aligningSingle01[2] + 1;
							/*if (scaffoldIndex == scaffoldIndexNum) {
								fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R0\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);
							}*/
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,R\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[6]);

						}
						//cout << endl;

					}

				}
				gapIndex = -1;
				gapInScaffoldIndex = -1;
				scaffoldIndex = -1;
				leftContigDiff = -1;
				rightContigDiff = -1;
				leftContigOverlap = -1;
				rightContigOverlap = -1;
				if (m+1 == num - 1) {//计算最后一个contig右边gap的对应序列
					
					int t = m + 1;
					gapIndex = GetGapIndexFromContigIndexOne(scaffoldSetHead, aligningSingle01[t* dataNum], 1, scaffoldIndex, gapInScaffoldIndex);
					//cout << "LscaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << endl;
					if (gapIndex != -1) {
						int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
						long int contigLength = contigSetHead->contigSet[aligningSingle01[t * dataNum]].contigLength;
						//cout << "LreadLength=" << readLength << "       contigLength=" << contigLength << "   Diff=" << contigLength - aligningSingle01[t * dataNum + 4] - 1 << endl;
						if (aligningSingle01[m * dataNum + 6]== aligningSingle01[t * dataNum + 6] && aligningSingle01[t * dataNum+6] == 1 && contigLength - aligningSingle01[t * dataNum+4] - 1 < readLength - aligningSingle01[t * dataNum+2] - 1 && contigLength - aligningSingle01[t * dataNum+4] - 1 < 600 && readLength - aligningSingle01[t * dataNum + 2] - 1>100) {
							//cout << "L1************************************************************************contigIndex=" << aligningSingle01[t* dataNum] << endl;
							int start = 0;
							int end = 0;
							leftContigDiff = contigLength - aligningSingle01[t * dataNum+4] - 1;
							leftContigOverlap = aligningSingle01[t * dataNum+5];
							if (readLength - aligningSingle01[t * dataNum+2] - 1 > distance + diff) {
								end = aligningSingle01[t * dataNum+2] + distance + diff;

							}
							else {
								end = readLength - 1;
							}
							start = aligningSingle01[t * dataNum+2] + 1;
							/*if (scaffoldIndex == scaffoldIndexNum) {
								fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L1\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[t * dataNum + 6]);
							}*/
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[t * dataNum+6]);
						}
						if (aligningSingle01[m * dataNum + 6] == aligningSingle01[t * dataNum + 6] && aligningSingle01[t * dataNum+6] == 0 && contigLength - aligningSingle01[t * dataNum+4] - 1 < aligningSingle01[t * dataNum+1] && contigLength - aligningSingle01[t * dataNum+4] - 1 < 600 && aligningSingle01[t * dataNum + 1]>100) {
							//cout << "L2************************************************************************contigIndex=" << aligningSingle01[t* dataNum] << endl;
							int start = 0;
							int end = 0;
							leftContigDiff = contigLength - aligningSingle01[t * dataNum+4] - 1;
							leftContigOverlap = aligningSingle01[t * dataNum+5];
							if (aligningSingle01[t * dataNum+1] > distance + diff) {
								start = aligningSingle01[t * dataNum+1] - distance - diff;
							}
							else {
								start = 0;
							}
							end = aligningSingle01[t * dataNum+1] - 1;
							/*if (scaffoldIndex == scaffoldIndexNum) {
								fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L1\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[t * dataNum + 6]);
							}*/
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,L\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[t * dataNum+6]);
						}
						//cout << endl;
					}

				}

				gapIndex = -1;
				gapInScaffoldIndex = -1;
				scaffoldIndex = -1;
				leftContigDiff = -1;
				rightContigDiff = -1;
				leftContigOverlap = -1;
				rightContigOverlap = -1;
				
				//方向为正的跨越读数-正向
				//cout << contigSetHead->contigSet[aligningSingle01[m * dataNum]].contigLength << "----------" << contigSetHead->contigSet[aligningSingle01[n * dataNum]].contigLength << endl;
				if (abs(aligningSingle01[m * dataNum] - aligningSingle01[n * dataNum]) == 1 && aligningSingle01[m * dataNum + 6] == aligningSingle01[n * dataNum + 6]
					&& aligningSingle01[m * dataNum + 6] == 1 && aligningSingle01[m * dataNum + 1] < aligningSingle01[n * dataNum + 1]) {
					long int contigLength = contigSetHead->contigSet[aligningSingle01[m * dataNum]].contigLength;

					//cout << "contig1=" << aligningSingle01[m * dataNum] << "--------" << "contig2=" << aligningSingle01[n * dataNum] << endl;
					//cout << contigSetHead->contigSet[aligningSingle01[m * dataNum]].contigLength << "----------" << contigSetHead->contigSet[aligningSingle01[n * dataNum]].contigLength << endl;
					//cout << "leftDiff=" << contigLength - aligningSingle01[m * dataNum + 4] - 1 << "          rightDiff=" << aligningSingle01[n * dataNum + 3] << endl;

					if (contigLength - aligningSingle01[m * dataNum + 4] - 1 < 8000 && aligningSingle01[n * dataNum + 3] < 8000) {
						gapIndex = GetGapIndexFromContigIndex(scaffoldSetHead, aligningSingle01[m * dataNum], aligningSingle01[n * dataNum], scaffoldIndex, gapInScaffoldIndex);
						/*if (scaffoldIndex == scaffoldIndexNum) {
							cout << "scaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << "-----------" << "gapLength=" << aligningSingle01[n * dataNum + 1] - aligningSingle01[m * dataNum + 2]+1 << "------initGapLength=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance << endl;
						}*/
					}
					if (gapIndex != -1) {
						int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
						if ((abs(distance - (aligningSingle01[n * dataNum + 1] - aligningSingle01[m * dataNum + 2]) < distance * 0.5) || abs(distance - (aligningSingle01[n * dataNum + 1] - aligningSingle01[m * dataNum + 2])) < 10000) && aligningSingle01[n * dataNum + 1] > aligningSingle01[m * dataNum + 2]) {
							/*if (scaffoldIndex == scaffoldIndexNum) {
								cout << "scaffoldIndex=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << "-----------" << "gapLength=" << aligningSingle01[n * dataNum + 1] - aligningSingle01[m * dataNum + 2] + 1 << "------initGapLength=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance << endl;
							}*/
							int start = -1;
							int end = -1;
							leftContigDiff = contigLength - aligningSingle01[m * dataNum + 4] - 1;
							rightContigDiff = aligningSingle01[n * dataNum + 3];
							
							leftContigOverlap = aligningSingle01[m * dataNum + 5];
							rightContigOverlap= aligningSingle01[n * dataNum + 5];

							start = aligningSingle01[m * dataNum + 2] + 1- extendReadLength;
							end = aligningSingle01[n * dataNum + 1] - 1+ extendReadLength;
							/*if (scaffoldIndex == scaffoldIndexNum) {
								fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,M\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end,leftContigDiff,rightContigDiff ,leftContigOverlap,rightContigOverlap,weigth, aligningSingle01[m * dataNum + 6]);
							}*/
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,M\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end,leftContigDiff,rightContigDiff ,leftContigOverlap,rightContigOverlap,weigth, aligningSingle01[m * dataNum + 6]);

						}
					}
				}

				//方向为负的跨越读数
				if (abs(aligningSingle01[m * dataNum] - aligningSingle01[n * dataNum]) == 1 && aligningSingle01[m * dataNum + 6] == aligningSingle01[n * dataNum + 6]
					&& aligningSingle01[m * dataNum + 6] == 0 && aligningSingle01[m * dataNum + 1] > aligningSingle01[n * dataNum + 1]) {
					long int contigLength = contigSetHead->contigSet[aligningSingle01[m * dataNum]].contigLength;

					//cout << "contig3=" << aligningSingle01[m * dataNum] << "--------" << "contig4=" << aligningSingle01[n * dataNum] << endl;
					//cout << contigSetHead->contigSet[aligningSingle01[m * dataNum]].contigLength << "----------" << contigSetHead->contigSet[aligningSingle01[n * dataNum]].contigLength << endl;
					
					//cout << "leftDiff=" << contigLength - aligningSingle01[m * dataNum + 4] - 1 << "          rightDiff=" << aligningSingle01[n * dataNum + 3] << endl;

					if (contigLength - aligningSingle01[m * dataNum + 4] - 1 < 500 && aligningSingle01[n * dataNum + 3] < 500) {
						gapIndex = GetGapIndexFromContigIndex(scaffoldSetHead, aligningSingle01[m * dataNum], aligningSingle01[n * dataNum], scaffoldIndex, gapInScaffoldIndex);
						/*if (scaffoldIndex == scaffoldIndexNum) {
							cout << "scaffoldIndex01=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << "-----------" << "gapLength=" << aligningSingle01[m * dataNum + 1] - aligningSingle01[n * dataNum + 2]+1 << "------initGapLength=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance << endl;
							for (long int h = 0; h < num; h++) {
								cout << aligningSingle01[h * dataNum] << "---"<<aligningSingle01[h * dataNum + 1] << "---"<<aligningSingle01[h * dataNum + 2] << "---"<<aligningSingle01[h * dataNum + 3] << "---"<<aligningSingle01[h * dataNum + 4] << "---"<<aligningSingle01[h * dataNum + 5] << "---"<<aligningSingle01[h * dataNum + 6] << endl;

							}
						}*/
						
					}
					if (gapIndex != -1) {
						int distance = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance;
						if ((abs(distance - (aligningSingle01[m * dataNum + 1] - aligningSingle01[n * dataNum + 2]) < distance * 0.5) || abs(distance - (aligningSingle01[m * dataNum + 1] - aligningSingle01[n * dataNum + 2])) < 10000) && aligningSingle01[m * dataNum + 1] > aligningSingle01[n * dataNum + 2]) {
							/*if (scaffoldIndex == scaffoldIndexNum) {
								cout << "scaffoldIndex02=" << scaffoldIndex << "**********" << "gapInScaffoldIndex=" << gapInScaffoldIndex << "-----------" << "gapLength=" << aligningSingle01[m * dataNum + 1] - aligningSingle01[n * dataNum + 2] + 1 << "------initGapLength=" << scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance << endl;
							
							}*/
							//cout << contigSetHead->contigSet[aligningSingle01[m * dataNum]].contigLength << "----------" << contigSetHead->contigSet[aligningSingle01[n * dataNum]].contigLength << endl;
							int start = -1;
							int end = -1;

							leftContigDiff = contigLength - aligningSingle01[m * dataNum + 4] - 1;
							rightContigDiff = aligningSingle01[n * dataNum + 3];

							leftContigOverlap = aligningSingle01[m * dataNum + 5];
							rightContigOverlap = aligningSingle01[n * dataNum + 5];

							start = aligningSingle01[n * dataNum + 2] + 1- extendReadLength;
							end = aligningSingle01[m * dataNum + 1] - 1+ extendReadLength;
							/*if (scaffoldIndex == scaffoldIndexNum) {
								fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,M\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap, weigth, aligningSingle01[m * dataNum + 6]);
							}*/
							fprintf(hifiFileFP, "%ld,%ld,%ld,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,M\n", scaffoldIndex, gapInScaffoldIndex, readIndex, start, end, leftContigDiff, rightContigDiff, leftContigOverlap, rightContigOverlap,weigth, aligningSingle01[m * dataNum + 6]);
						}
					}
				}


			}
			num = 0;
		}
	}
	
	fflush(hifiFileFP);
	fclose(fp);
	
}

//表示contigIndex的左边gap索引，或者右边gap索引
long int GetGapIndexFromContigIndexOne(ScaffoldSetHead* scaffoldSetHead, long int contigIndex, bool right, long int& scaffoldIndex, long int& gapIndex) {
	if (contigIndex < 0) {
		return -1;
	}
	long int tempContigCount = 0;
	long int tempGapIndex = 0;
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		if (tempContigCount + scaffoldSetHead->scaffoldSet[i].contigCount <= contigIndex) {
			tempContigCount = tempContigCount + scaffoldSetHead->scaffoldSet[i].contigCount;
			tempGapIndex = tempGapIndex + scaffoldSetHead->scaffoldSet[i].gapCount;
			continue;
		}
		int temp = 0;
		if (right == true) {//right非0
			temp = contigIndex - tempContigCount;
		}
		else {//right为0
			temp = contigIndex - tempContigCount - 1;
		}
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {
			temp++;
		}
		if (temp < 0 || temp >= scaffoldSetHead->scaffoldSet[i].gapCount) {
			scaffoldIndex = -1;
			gapIndex = -1;
			return -1;
		}
		scaffoldIndex = i;
		gapIndex = temp;
		tempGapIndex = tempGapIndex + gapIndex;
		return tempGapIndex;
	}
	scaffoldIndex = -1;
	gapIndex = -1;
	return -1;
}

//通过两个相邻的contig索引确定相应的gap
long int GetGapIndexFromContigIndex(ScaffoldSetHead * scaffoldSetHead, long int leftContigIndex, long int rightContigIndex, long int& scaffoldIndex, long int& gapIndex) {
	if (abs(leftContigIndex - rightContigIndex) != 1 || rightContigIndex < 0 || leftContigIndex < 0) {
		scaffoldIndex = -1;
		gapIndex = -1;
		return -1;
	}
	if (leftContigIndex > rightContigIndex) {
		long int temp = leftContigIndex;
		leftContigIndex = rightContigIndex;
		rightContigIndex = temp;
	}
	long int tempContigCount = 0;
	long int tempGapIndex = 0;
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		if (tempContigCount + scaffoldSetHead->scaffoldSet[i].contigCount <= leftContigIndex) {
			tempContigCount = tempContigCount + scaffoldSetHead->scaffoldSet[i].contigCount;
			tempGapIndex = tempGapIndex + scaffoldSetHead->scaffoldSet[i].gapCount;
			continue;
		}
		if (tempContigCount + scaffoldSetHead->scaffoldSet[i].contigCount == rightContigIndex) {
			break;
		}
		int temp = leftContigIndex - tempContigCount;
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {
			temp++;
		}
		//cout << "temp=" << temp << endl;
		scaffoldIndex = i;
		gapIndex = temp;
		tempGapIndex = tempGapIndex + gapIndex;
		return tempGapIndex;
	}
	scaffoldIndex = -1;
	gapIndex = -1;
	return -1;
}

//通过gap计算gap两侧的contig
long int GetContigIndexFromGapIndex(ScaffoldSetHead* scaffoldSetHead, long int scaffoldIndex,long int gapIndex,long int &gapLeftContigindex,long int &gapRightContigIndex) {
	long int tempContigCount = 0;
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		if (i < scaffoldIndex) {
			tempContigCount = tempContigCount + scaffoldSetHead->scaffoldSet[i].contigCount;
			continue;
		}
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {
			if (gapIndex == 0) {
				gapLeftContigindex = -1;
				gapRightContigIndex = tempContigCount + gapIndex + 1;
			}
			else {
				gapLeftContigindex = tempContigCount + gapIndex;
				gapRightContigIndex = tempContigCount + gapIndex + 1;
			}
			
		}else {
			gapLeftContigindex = tempContigCount + gapIndex + 1;
			gapRightContigIndex = tempContigCount + gapIndex + 2;
		}
		return 1;
	}
	return -1;
}

//获取对应gap的三种类型的序列
void GetHiFiInGap(ContigSetHead* contigSetHead, ScaffoldSetHead* scaffoldSetHead, char* hifiFile, char* positionFile, char* line, int maxSize, char* outputAddress) {
	long int diff =200;
	long int errorDiff = 500;

	FILE* fp;
	if ((fp = fopen(positionFile, "r")) == NULL) {
		printf("%s, does not exist!", positionFile);
		exit(0);
	}

	FILE* fp1;
	if ((fp1 = fopen(hifiFile, "r")) == NULL) {
		printf("%s, does not exist!", hifiFile);
		exit(0);
	}
	long int hifiCount = 0;

	while ((fgets(line, maxSize, fp1)) != NULL) {
		if (line[0] == '>') {
			hifiCount++;
		}
	}
	fclose(fp1);
	//cout << hifiCount << endl;

	bool* hifiGapTrue = (bool*)malloc(sizeof(bool) * hifiCount);//用于标记使用哪些读数
	for (long int i = 0; i < hifiCount; i++) {
		hifiGapTrue[i] = false;
	}
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		scaffoldSetHead->scaffoldSet[i].hifiGapSet = (HifiGapSet *)malloc(sizeof(HifiGapSet) * scaffoldSetHead->scaffoldSet[i].gapCount);//只有有gap的scaffold才会malloc
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {

			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet = NULL;//相应读数碱基
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet = NULL;
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet = NULL;

			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount = 0;//读数数量
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount = 0;
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount = 0;

			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftTooLengthCount = 0;
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightTooLengthCount = 0;

			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence = NULL;//一致序列
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftConsensusSequence = NULL;
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightConsensusSequence = NULL;
		}
	}
	const char* split = ",";
	char* p;
	long int scaffoldIndex = -1;
	long int gapInScaffoldIndex = -1;
	long int readIndex = -1;
	long int readStartPosition = -1;
	long int readEndPosition = -1;
	long int leftContigDiff = -1;
	long int rightContigDiff = -1;
	long int leftContigOverlap = 0;
	long int rightContigOverlap = 0;
	long int weigth = 0;
	bool orientation = 0;
	
	long int initialReadLength = 0;

	while ((fgets(line, maxSize, fp)) != NULL) {//统计每个gap对应的三类读数的数量，并标记使用的读数
		p = strtok(line, split);
		scaffoldIndex = atoi(p);

		p = strtok(NULL, split);
		gapInScaffoldIndex = atoi(p);

		p = strtok(NULL, split);
		readIndex = atoi(p);
		hifiGapTrue[readIndex] = true;//确定需要的HIFI读数

		p = strtok(NULL, split);
		readStartPosition = atoi(p);

		p = strtok(NULL, split);
		readEndPosition = atoi(p);

		p = strtok(NULL, split);
		leftContigDiff = atoi(p);

		p = strtok(NULL, split);
		rightContigDiff = atoi(p);

		p = strtok(NULL, split);
		leftContigOverlap = atoi(p);

		p = strtok(NULL, split);
		rightContigOverlap = atoi(p);

		p = strtok(NULL, split);
		weigth = atoi(p);

		p = strtok(NULL, split);
		orientation = atoi(p);

		p = strtok(NULL, split);

		//cout << scaffoldIndex << "    " << gapInScaffoldIndex << "    " << readIndex << "   " << readStartPosition << "  " << readEndPosition << "  " << leftContigDiff << "    " << rightContigDiff << "   " << leftContigOverlap << "   " << rightContigOverlap << "    " << weigth << "    " << orientation << "    " << p[0] << endl;

		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance <= 0) {
			continue;
		}

		if (p[0] == 'M') {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadCount++;
		}
		else if(p[0]=='L') {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadCount++;
		}
		else {
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadCount++;
		}

	}
	fclose(fp);


	ReadSetHead* readSetHead = GetReadSet(hifiFile,hifiCount, hifiGapTrue);//只获取需要的读数
	
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {//为存储gap对应序列创建空间
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {//j表示第几个gap
			
			/*cout << " i=" << i << "  j=" << j << "  readCount:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount << "---- " << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount << "---" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount
				<< "       estimatedGapDistance=" << scaffoldSetHead->scaffoldSet[i].gapSet[j].estimatedGapDistance << "      gapLength=" << scaffoldSetHead->scaffoldSet[i].gapSet[j].gapLength << endl;*/

			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount != 0) {
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet = (ReadToGap*)malloc(sizeof(ReadToGap) * scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount);
				
			}
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount != 0) {
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet = (ReadToGap*)malloc(sizeof(ReadToGap) * scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount);
			}
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount != 0) {
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet = (ReadToGap*)malloc(sizeof(ReadToGap) * scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount);
			}
		}
		//cout << endl;
	}

	if ((fp = fopen(positionFile, "r")) == NULL) {
		printf("%s, does not exist!", positionFile);
		exit(0);
	}

	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		//for (long int k = 0; k < scaffoldSetHead->scaffoldSet[i].gapCount; k++) {//输出每个gap对应的三类读数数量
		//	cout << "i=" << i << "   k=" << k << "-" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[k].spanReadCount << "-" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[k].leftReadCount << "-" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[k].rightReadCount << endl;
		//}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount = 0;
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount = 0;
			scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount = 0;
		}
	}
	long int contigIndex = -1;
	
	while ((fgets(line, maxSize, fp)) != NULL) {
	/*for(int k=0;k<10;k++){
		fgets(line, maxSize, fp);*/
		p = strtok(line, split);
		scaffoldIndex = atoi(p);

		p = strtok(NULL, split);
		gapInScaffoldIndex = atoi(p);

		p = strtok(NULL, split);
		readIndex = atoi(p);

		p = strtok(NULL, split);
		readStartPosition = atoi(p);

		p = strtok(NULL, split);
		readEndPosition = atoi(p);

		p = strtok(NULL, split);
		leftContigDiff = atoi(p);

		p = strtok(NULL, split);
		rightContigDiff = atoi(p);

		p = strtok(NULL, split);
		leftContigOverlap = atoi(p);

		p = strtok(NULL, split);
		rightContigOverlap = atoi(p);

		p = strtok(NULL, split);
		weigth = atoi(p);

		p = strtok(NULL, split);
		orientation = atoi(p);

		p = strtok(NULL, split);

		//cout << scaffoldIndex << "    " << gapInScaffoldIndex <<"    " <<readIndex<<"   "<< readStartPosition<<"  "<< readEndPosition <<"  "<< leftContigDiff <<"    "<< rightContigDiff<<"   "<< leftContigOverlap <<"   "<< rightContigOverlap << "    " << weigth << "    " << orientation << "    " << p[0] << endl;


		if (scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapInScaffoldIndex].estimatedGapDistance <= 0) {
			continue;
		}

		char* tempRead = (char*)malloc(sizeof(char) * (readEndPosition - readStartPosition + 2));
		initialReadLength = strlen(readSetHead->readSet[readIndex].read);//原始读数长度
		//cout << "initialReadLength=" << initialReadLength << endl;
		strncpy(tempRead, readSetHead->readSet[readIndex].read + readStartPosition, readEndPosition - readStartPosition + 1);
		tempRead[readEndPosition - readStartPosition + 1] = '\0';
		if (orientation != true) {
			tempRead = ReverseComplement(tempRead);
		}

		if (p[0] == 'M') {
			long int index = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadCount++;

			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadSet[index].read = tempRead;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadSet[index].readLength = readEndPosition - readStartPosition + 1;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadSet[index].leftContigDiff = leftContigDiff;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadSet[index].rightContigDiff = rightContigDiff;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadSet[index].leftContigOverlap = leftContigOverlap;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadSet[index].rightContigOverlap = rightContigOverlap;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].spanReadSet[index].weigth = weigth;
		}
		else if (p[0] == 'L') {
			long int index = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadCount++;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadSet[index].read = tempRead;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadSet[index].readLength = readEndPosition - readStartPosition + 1;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadSet[index].leftContigDiff=leftContigDiff;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadSet[index].rightContigDiff= rightContigDiff;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadSet[index].leftContigOverlap= leftContigOverlap;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadSet[index].rightContigOverlap= rightContigOverlap;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftReadSet[index].weigth= weigth;
			
			if (orientation == true) {
				if (readEndPosition< initialReadLength - 1) {
					if (initialReadLength -1- readEndPosition+diff > errorDiff) {
						scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftTooLengthCount += 1;
					}
				}
			}
			else {
				if (readStartPosition+diff> errorDiff) {
					scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].leftTooLengthCount += 1;
				}
			}
		}
		else {
			long int index = scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadCount++;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadSet[index].read = tempRead;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadSet[index].readLength = readEndPosition - readStartPosition + 1;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadSet[index].leftContigDiff = leftContigDiff;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadSet[index].rightContigDiff = rightContigDiff;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadSet[index].leftContigOverlap = leftContigOverlap;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadSet[index].rightContigOverlap = rightContigOverlap;
			scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightReadSet[index].weigth = weigth;

			if (orientation == true) {
				if (readStartPosition+diff > errorDiff) {
					scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightTooLengthCount += 1;
				}
			}
			else {
				if (initialReadLength-1-readEndPosition+diff> errorDiff) {
					scaffoldSetHead->scaffoldSet[scaffoldIndex].hifiGapSet[gapInScaffoldIndex].rightTooLengthCount += 1;
				}
			}
		}
	}

	fclose(fp);
	long int gapLeftContigindex = -1;
	long int gapRightContigIndex = -1;

	long int m = 0;
	long int n = 0;
	long int index01 = 0;
	char* read01;
	long int readLength01=0;
	long int weigth01 = 0;
	long int leftContigDiff01 = 0;
	long int rightContigDiff01 = 0;
	long int leftContigOverlap01 = 0;
	long int rightContigOverlap01 = 0;

	
	
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {//为spanread按长度排序
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount < 2) {
				continue;
			}
			for (m = 0; m < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount-1; m++) {
				readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].readLength;
				index01 = m;
				for (n = m + 1; n < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount; n++) {
					
					if (readLength01 > scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[n].readLength) {
						index01 = n;
						readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[n].readLength;
					}
				}
				if (index01 == m) {
					continue;
				}
				read01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].read;
				readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].readLength;
				weigth01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].weigth;
				leftContigDiff01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].leftContigDiff;
				rightContigDiff01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].rightContigDiff;
				leftContigOverlap01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].leftContigOverlap;
				rightContigOverlap01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].rightContigOverlap;

				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].read = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].read;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].readLength = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].readLength;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].weigth = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].weigth;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].leftContigDiff = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].leftContigDiff;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].rightContigDiff = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].rightContigDiff;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].leftContigOverlap = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].leftContigOverlap;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].rightContigOverlap = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].rightContigOverlap;

				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].read = read01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].readLength = readLength01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].weigth = weigth01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].leftContigDiff = leftContigDiff01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].rightContigDiff = rightContigDiff01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].leftContigOverlap = leftContigOverlap01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[index01].rightContigOverlap = rightContigOverlap01;
			}
		}
		
	}

	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {//leftread按长度排序
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount < 2) {
				continue;
			}
			for (m = 0; m < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount - 1; m++) {
				readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].readLength;
				index01 = m;
				for (n = m + 1; n < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount; n++) {

					if (readLength01 > scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[n].readLength) {
						index01 = n;
						readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[n].readLength;
					}
				}
				if (index01 == m) {
					continue;
				}
				read01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].read;
				readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].readLength;
				weigth01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].weigth;
				leftContigDiff01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].leftContigDiff;
				rightContigDiff01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].rightContigDiff;
				leftContigOverlap01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].leftContigOverlap;
				rightContigOverlap01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].rightContigOverlap;

				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].read = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].read;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].readLength = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].readLength;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].weigth = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].weigth;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].leftContigDiff = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].leftContigDiff;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].rightContigDiff = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].rightContigDiff;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].leftContigOverlap = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].leftContigOverlap;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[m].rightContigOverlap = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].rightContigOverlap;

				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].read = read01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].readLength = readLength01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].weigth = weigth01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].leftContigDiff = leftContigDiff01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].rightContigDiff = rightContigDiff01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].leftContigOverlap = leftContigOverlap01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadSet[index01].rightContigOverlap = rightContigOverlap01;
			}
		}

	}
	
	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {//rightRead按长度排序
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount < 2) {
				continue;
			}
			for (m = 0; m < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount - 1; m++) {
				readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].readLength;
				index01 = m;
				for (n = m + 1; n < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount; n++) {

					if (readLength01 > scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[n].readLength) {
						index01 = n;
						readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[n].readLength;
					}
				}
				if (index01 == m) {
					continue;
				}
				read01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].read;
				readLength01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].readLength;
				weigth01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].weigth;
				leftContigDiff01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].leftContigDiff;
				rightContigDiff01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].rightContigDiff;
				leftContigOverlap01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].leftContigOverlap;
				rightContigOverlap01 = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].rightContigOverlap;

				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].read = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].read;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].readLength = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].readLength;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].weigth = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].weigth;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].leftContigDiff = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].leftContigDiff;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].rightContigDiff = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].rightContigDiff;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].leftContigOverlap = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].leftContigOverlap;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[m].rightContigOverlap = scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].rightContigOverlap;

				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].read = read01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].readLength = readLength01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].weigth = weigth01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].leftContigDiff = leftContigDiff01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].rightContigDiff = rightContigDiff01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].leftContigOverlap = leftContigOverlap01;
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadSet[index01].rightContigOverlap = rightContigOverlap01;
			}
		}

	}

	//for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {//为存储gap对应序列创建空间
	//	if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
	//		continue;
	//	}
	//	if (i != 6649) {
	//		continue;
	//	}
	//	for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
	//		if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount < 2) {
	//			continue;
	//		}
	//		if (j != 12) {
	//			continue;
	//		}
	//		for (m = 0; m < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount; m++) {
	//			cout << "i=" << i << "   j=" << j << "   m=" << m << endl;
	//			cout << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].weigth << "----Diff:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].leftContigDiff
	//				<< "---" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].rightContigDiff << "---Overlap:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].leftContigOverlap << "----" 
	//				<< scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].rightContigOverlap <<"*********************************readLength="<< scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[m].readLength<< endl;

	//		}
	//		cout << endl;
	//	}
	//	
	//}
	
	/*for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		if (i > 6644) {
			continue;
		}

		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
			cout << ">" << i << "--" << j << "----readCount:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount << "-------" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount << "-------" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount << endl;
			
			//for (long int k = 0; k < scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount; k++) {//输出gap对应的读数
			//	cout << ">" << i <<"-------" <<j<<"    " << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[k].leftContigDiff << "----" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[k].leftContigDiff
			//		<<"-----"<<scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[k].leftContigOverlap << "----" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[k].leftContigOverlap 
			//		<<"-------"<< scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[k].readLength << endl;
			//	
			//	
			//	
			//	cout << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[k].read << endl;
			//	//cout << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadSet[k].readLength << endl;
			//}
			//cout <<"i="<<i<<"    j="<<j<<"   contigLength="<< scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength << endl;



			

		}
		cout << endl;

	}*/
	

	long int a = 6658;//scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadCount
	long int b = 4;
	//for (long int k = 0; k < scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadCount; k++) {//输出gap对应的读数
	//	cout << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadCount << endl;
	//	cout << ">M" << k << "--" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadSet[k].weigth << "-------" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadSet[k].readLength <<
	//		"-----gap:" << scaffoldSetHead->scaffoldSet[a].gapSet[b].estimatedGapDistance << "-----" << scaffoldSetHead->scaffoldSet[a].gapSet[b].gapLength<<
	//		"-----overlap:"<< scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadSet[k].leftContigOverlap<<"------"<< scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadSet[k].rightContigOverlap<<
	//		"-----diff:"<< scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadSet[k].leftContigDiff<<"-----" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadSet[k].rightContigDiff << endl;
	//	cout << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].spanReadSet[k].read << endl;
	//}

	//for (long int k = 0; k < scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadCount; k++) {//输出gap对应的读数
	//	cout << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadCount << endl;
	//	cout << ">L" << k << "--" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadSet[k].weigth << "-------" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadSet[k].readLength <<
	//		"-----gap:" << scaffoldSetHead->scaffoldSet[a].gapSet[b].estimatedGapDistance << "-----" << scaffoldSetHead->scaffoldSet[a].gapSet[b].gapLength <<
	//		"-----overlap:" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadSet[k].leftContigOverlap << "------" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadSet[k].rightContigOverlap <<
	//		"-----diff:" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadSet[k].leftContigDiff << "-----" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadSet[k].rightContigDiff << endl;
	//	cout << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].leftReadSet[k].read << endl;
	//}

	//for (long int k = 0; k < scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadCount; k++) {//输出gap对应的读数
	//	cout << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadCount << endl;
	//	cout << ">R" << k << "--" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadSet[k].weigth << "-------" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadSet[k].readLength <<
	//		"-----gap:" << scaffoldSetHead->scaffoldSet[a].gapSet[b].estimatedGapDistance << "-----" << scaffoldSetHead->scaffoldSet[a].gapSet[b].gapLength <<
	//		"-----overlap:" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadSet[k].leftContigOverlap << "------" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadSet[k].rightContigOverlap <<
	//		"-----diff:" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadSet[k].leftContigDiff << "-----" << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadSet[k].rightContigDiff << endl;
	//	cout << scaffoldSetHead->scaffoldSet[a].hifiGapSet[b].rightReadSet[k].read << endl;
	//}

	

}

















#endif