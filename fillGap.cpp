#ifndef FILLGAP_CPP_INCLUDED 
#define FILLGAP_CPP_INCLUDED 

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
#include"raconConsensus.h"
#include"fillGap.h"

using namespace std;

void FillingGap(ScaffoldSetHead* scaffoldSetHead, char* resultOutPutDirectory) {
	
    long int gapNoNullCount = 0;
    long int combinSpanCount = 0;
	long int gapIsNullCount = 0;
	long int totalGapCount = 0;

	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			continue;
		}
		/*if (i != 6637) {
			continue;
		}*/

		totalGapCount = totalGapCount + scaffoldSetHead->scaffoldSet[i].gapCount;
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
			/*if (j != 1){
				continue;
			}*/
			cout << endl;
			cout << "********************************************************************A gap start***********************************************************************************" << endl;
			cout << "start-filling:" << i << "--" << j << "--"<< scaffoldSetHead->scaffoldSet[i].gapSet[j].gapLength << endl;
			cout << "readCount:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanReadCount << "--" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftReadCount << "--" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightReadCount << endl;
			cout<<"extendContigNum:"<< scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftTooLengthCount << "---" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightTooLengthCount << endl;
			
			bool result=false;

			if (result == false) {
				result = GetRaconConsensusSequenceForOneGap(scaffoldSetHead, i, j,resultOutPutDirectory);
				//result = GetMafftConsensusSequenceForOneGap(scaffoldSetHead, i, j, resultOutPutDirectory);
			}
			
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence != NULL) {
				gapNoNullCount++;
				cout << "finalConsensus:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence << endl;
			}
			
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence == NULL) {
				gapIsNullCount++;
				
				scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence = CombineLeftAndRightConsensusSequence(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftConsensusSequence,scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightConsensusSequence, scaffoldSetHead->scaffoldSet[i].gapSet[j].gapLength);
				if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence != NULL) {
					combinSpanCount++;
					cout << "combineConsensus:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence << endl;
				}
			}

			//TransferToLower(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence);
			//cout << endl;
			cout << "********************************************************************A gap End***********************************************************************************" << endl;
			cout << endl;
		}
		
	}

	//for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
	//	if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
	//		continue;
	//	}
	//	
	//	for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].gapCount; j++) {
	//		if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence != NULL) {
	//			cout<<"S>"<<i<<"     j="<<j<< "   readLength=" << strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence) <<"    gapLength=" << scaffoldSetHead->scaffoldSet[i].gapSet[j].gapLength << "  estimatedGapDistance=" << scaffoldSetHead->scaffoldSet[i].gapSet[j].estimatedGapDistance << endl;
	//			cout << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].spanConsensusSequence << endl;
	//		}
	//		/*if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightConsensusSequence != NULL) {
	//			cout << "L>" << i << "     j=" << j << endl;
	//			cout << "finalRightConsensus:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].rightConsensusSequence << endl;
	//		}
	//		if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftConsensusSequence != NULL) {
	//			cout << "R>" << i << "     j=" << j << endl;
	//			cout << "finalLeftConsensus:" << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j].leftConsensusSequence << endl;
	//		}*/
	//	}
	//}





	char fillingResultFile[1000];
	sprintf(fillingResultFile, "%s/gapfilling-scaffoldset.fa", resultOutPutDirectory);
	OutputFillingResultFromScaffoldSetHead1(scaffoldSetHead, fillingResultFile);
	cout << "totalGapCount="<< totalGapCount<<"-------gapNoNullCount=" << gapNoNullCount << "--------" << "gapIsNullCount=" << gapIsNullCount <<"----------"<<"combinSpanCount=" << combinSpanCount<<endl;
	
	cout << "FillingGap.cpp" << endl;
	
}


void OutputFillingResultFromScaffoldSetHead1(ScaffoldSetHead* scaffoldSetHead, char* fillingResultFile) {
	double indentityRatio = 0.98;

	FILE* fp;
	if ((fp = fopen(fillingResultFile, "w")) == NULL) {
		printf("%s, does not exist!", fillingResultFile);
		exit(0);
	}
	long int gapIndex = 0;
	long int length = 1000000;
	char* contig = (char*)malloc(sizeof(char) * length);
	int flag1 = 0;

	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {
			fprintf(fp, ">%ld\n%s\n", i + flag1, scaffoldSetHead->scaffoldSet[i].scaffold);
			continue;
		}
		fprintf(fp, ">%ld\n", i + flag1);
		gapIndex = 0;
		int flag = 0;
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {//
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence == NULL) {
				long int t = 0;
				for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[0].gapLength; t++) {
					contig[t] = 'N';
				}
				contig[t] = '\0';
			}
			else {
				strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence);
				contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence)] = '\0';
			}

			fprintf(fp, "%s", contig);
			flag = 1;
			gapIndex++;//
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {//
			cout << "gapFilling:" << i << "-" << j << endl;
			if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength < length - 1) {//
				strncpy(contig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				contig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s", contig);
			}
			else {//

				char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength + 1);
				strncpy(tempContig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				tempContig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s", tempContig);
			}

			if (gapIndex < scaffoldSetHead->scaffoldSet[i].gapCount) {//

				if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanReadCount == 0) {
					/*
					if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength>4000 && scaffoldSetHead->scaffoldSet[i].contigPostion[j+1].contigLength>4000 &&
						scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftTooLengthCount >= scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftReadCount * 0.5 && scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightTooLengthCount >= scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightReadCount * 0.5
						&& (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftReadCount >0 && scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightReadCount>0)) {
						
						
						cout << "*******************************************************************breakPoint" << endl;
						if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftConsensusSequence != NULL) {
							if (strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftConsensusSequence) < length - 1) {
								strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftConsensusSequence);
								contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftConsensusSequence)] = '\0';
								fprintf(fp, "%s", contig);

							}
							else {
								char* tempContig = (char*)malloc(sizeof(char) * strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftConsensusSequence) + 1);
								strcpy(tempContig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftConsensusSequence);
								tempContig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].leftConsensusSequence)] = '\0';
								fprintf(fp, "%s", tempContig);

							}
						}
						
						fprintf(fp, "\n>%ld\n", i + flag1);
						flag1 = flag1 + 1;
						if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightConsensusSequence != NULL) {
							if (strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightConsensusSequence) < length - 1) {
								strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightConsensusSequence);
								contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightConsensusSequence)] = '\0';
								fprintf(fp, "%s", contig);

							}
							else {
								char* tempContig = (char*)malloc(sizeof(char) * strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightConsensusSequence) + 1);
								strcpy(tempContig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightConsensusSequence);
								tempContig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].rightConsensusSequence)] = '\0';
								fprintf(fp, "%s", tempContig);

							}
						}
						
					}
					else {
						if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence != NULL) {
							if (strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) < length - 1) {
								strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
								contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
								fprintf(fp, "%s", contig);

							}
							else {
								char* tempContig = (char*)malloc(sizeof(char) * strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) + 1);
								strcpy(tempContig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
								tempContig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
								fprintf(fp, "%s", tempContig);

							}

						}
						else {
							long int t = 0;
							for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength; t++) {
								contig[t] = 'N';
							}
							contig[t] = '\0';
							fprintf(fp, "%s", contig);

						}

					}
					*/
					if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence != NULL) {
						if (strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) < length - 1) {
							strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
							contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
							fprintf(fp, "%s", contig);

						}
						else {
							char* tempContig = (char*)malloc(sizeof(char) * strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) + 1);
							strcpy(tempContig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
							tempContig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
							fprintf(fp, "%s", tempContig);

						}

					}
					else {
						long int t = 0;
						for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength; t++) {
							contig[t] = 'N';
						}
						contig[t] = '\0';
						fprintf(fp, "%s", contig);

					}
				}
				else {
					if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence != NULL) {
						if (strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) < length - 1) {
							strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
							contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
							fprintf(fp, "%s", contig);

						}
						else {
							char* tempContig = (char*)malloc(sizeof(char) * strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) + 1);
							strcpy(tempContig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
							tempContig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
							fprintf(fp, "%s", tempContig);

						}

					}
					else {
						long int t = 0;
						for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength; t++) {
							contig[t] = 'N';
						}
						contig[t] = '\0';
						fprintf(fp, "%s", contig);

					}
				}
				
				gapIndex++;

			}
		}
		fprintf(fp, "\n");

	}
	fflush(fp);
	fclose(fp);

}





//��ͨ������
/*
void OutputFillingResultFromScaffoldSetHead1(ScaffoldSetHead* scaffoldSetHead, char* fillingResultFile) {
	double indentityRatio = 0.98;

	FILE* fp;
	if ((fp = fopen(fillingResultFile, "w")) == NULL) {
		printf("%s, does not exist!", fillingResultFile);
		exit(0);
	}
	long int gapIndex = 0;
	long int length = 1000000;
	char* contig = (char*)malloc(sizeof(char) * length);
	int flag1 = 0;//��ʾ��Ͽ�scaffold���ӵ�scaffold

	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {//�����scaffoldû��gapֱ�����
			fprintf(fp, ">%ld\n%s\n", i + flag1, scaffoldSetHead->scaffoldSet[i].scaffold);
			continue;
		}
		fprintf(fp, ">%ld\n", i + flag1);
		gapIndex = 0;
		int flag = 0;
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {//���scaffold�տ�ʼ����gap�����gap��Ӧһ�����з�������N
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence == NULL) {
				long int t = 0;
				for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[0].gapLength; t++) {
					contig[t] = 'N';
				}
				contig[t] = '\0';
			}
			else {
				strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence);
				contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence)] = '\0';
			}

			fprintf(fp, "%s", contig);
			flag = 1;
			gapIndex++;//��ʾ�ڼ���gap
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {//���scaffold�е�j��contig
			cout << "gapFilling:" << i << "-" << j << endl;
			if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength < length - 1) {//�����contig����С��length����ȡ��contig
				strncpy(contig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				contig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s", contig);
			}
			else {//�����contig���ȴ���length

				char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength + 1);
				strncpy(tempContig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				tempContig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s", tempContig);
			}

			if (gapIndex < scaffoldSetHead->scaffoldSet[i].gapCount) {//�������
				if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence != NULL) {
					if (strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) < length - 1) {
						strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
						contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
						fprintf(fp, "%s", contig);

					}
					else {
						char* tempContig = (char*)malloc(sizeof(char) * strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) + 1);
						strcpy(tempContig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
						tempContig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
						fprintf(fp, "%s", tempContig);

					}

				}
				else {
					long int t = 0;
					for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength; t++) {
						contig[t] = 'N';
					}
					contig[t] = '\0';
					fprintf(fp, "%s", contig);

				}
				gapIndex++;

			}
		}
		fprintf(fp, "\n");

	}
	fflush(fp);
	fclose(fp);

}*/


/*
void OutputFillingResultFromScaffoldSetHead1(ScaffoldSetHead* scaffoldSetHead, char* fillingResultFile) {
	double indentityRatio = 0.98;
	long int scaffoldCount = 6658;
	long int gapCount = 4;

	FILE* fp;
	if ((fp = fopen(fillingResultFile, "w")) == NULL) {
		printf("%s, does not exist!", fillingResultFile);
		exit(0);
	}
	long int gapIndex = 0;
	long int length = 1000000;
	char* contig = (char*)malloc(sizeof(char) * length);
	int flag1 = 0;//��ʾ��Ͽ�scaffold���ӵ�scaffold

	for (long int i = 0; i < scaffoldSetHead->scaffoldCount; i++) {
		if (scaffoldSetHead->scaffoldSet[i].gapCount < 1) {//�����scaffoldû��gapֱ�����
			fprintf(fp, ">%ld\n%s\n", i + flag1, scaffoldSetHead->scaffoldSet[i].scaffold);
			continue;
		}
		fprintf(fp, ">%ld\n", i + flag1);
		gapIndex = 0;
		int flag = 0;
		if (scaffoldSetHead->scaffoldSet[i].gapSet[0].gapStartCoordinate == 0) {//���scaffold�տ�ʼ����gap�����gap��Ӧһ�����з�������N
			if (scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence == NULL) {
				long int t = 0;
				for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[0].gapLength; t++) {
					contig[t] = 'N';
				}
				contig[t] = '\0';
			}
			else {
				strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence);
				contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[0].spanConsensusSequence)] = '\0';
			}

			fprintf(fp, "%s", contig);
			flag = 1;
			gapIndex++;//��ʾ�ڼ���gap
		}
		for (long int j = 0; j < scaffoldSetHead->scaffoldSet[i].contigCount; j++) {//���scaffold�е�j��contig
			cout << "gapFilling:" << i << "-" << j << endl;
			if (scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength < length - 1) {//�����contig����С��length����ȡ��contig
				strncpy(contig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				contig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s", contig);
			}
			else {//�����contig���ȴ���length

				char* tempContig = (char*)malloc(sizeof(char) * scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength + 1);
				strncpy(tempContig, scaffoldSetHead->scaffoldSet[i].scaffold + scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigStartCoordinate, scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength);
				tempContig[scaffoldSetHead->scaffoldSet[i].contigPostion[j].contigLength] = '\0';
				fprintf(fp, "%s", tempContig);
			}

			if (gapIndex < scaffoldSetHead->scaffoldSet[i].gapCount) {//�������
				if (i== scaffoldCount && j + flag==gapCount) {
					cout << scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence << endl;
					if (strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) < length - 1) {
						strcpy(contig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
						contig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
						fprintf(fp, "%s", contig);

					}
					else {
						char* tempContig = (char*)malloc(sizeof(char) * strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence) + 1);
						strcpy(tempContig, scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence);
						tempContig[strlen(scaffoldSetHead->scaffoldSet[i].hifiGapSet[j + flag].spanConsensusSequence)] = '\0';
						fprintf(fp, "%s", tempContig);

					}

				}
				else {
					long int t = 0;
					cout << "gapLength=" << scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength << endl;
					for (t = 0; t < scaffoldSetHead->scaffoldSet[i].gapSet[j + flag].gapLength; t++) {
						contig[t] = 'N';
					}
					contig[t] = '\0';
					fprintf(fp, "%s", contig);

				}
				gapIndex++;

			}
		}
		fprintf(fp, "\n");

	}
	fflush(fp);
	fclose(fp);

}
*/



//�ϲ���������

char* CombineLeftAndRightConsensusSequence(char* leftConsensusSequence, char* rightConsensusSequence, long int gapDistance) {
	long int gapIndex = 0;
	long int length = 0;
	long int length1 = 0;
	char* contig = NULL;
	long int dis = 5;
	if (leftConsensusSequence == NULL && rightConsensusSequence == NULL) {
		cout << "gap is NULL" << endl;
	}

	if (leftConsensusSequence != NULL && rightConsensusSequence == NULL) {
		length = strlen(leftConsensusSequence);
		if (length <= gapDistance) {
			contig = (char*)malloc(sizeof(char) * (gapDistance + 1));
			strncpy(contig, leftConsensusSequence, length);
			for (long int i = length; i < gapDistance; i++) {
				contig[i] = 'N';
			}
			contig[gapDistance] = '\0';
		}
		else {
			contig = (char*)malloc(sizeof(char) * (gapDistance + 1));
			strncpy(contig, leftConsensusSequence, gapDistance);
			contig[gapDistance] = '\0';
		}

	}

	if (leftConsensusSequence == NULL && rightConsensusSequence != NULL) {
		length = strlen(rightConsensusSequence);
		if (length <= gapDistance) {
			contig = (char*)malloc(sizeof(char) * (gapDistance + 1));
			strncpy(contig + gapDistance - length, rightConsensusSequence, length);
			for (long int i = 0; i < gapDistance - length; i++) {//�˴��޸Ľ�+1ȥ��
				contig[i] = 'N';
			}
			contig[gapDistance] = '\0';
		}
		else {
			contig = (char*)malloc(sizeof(char) * (gapDistance + 1));
			strncpy(contig, rightConsensusSequence + length - gapDistance, gapDistance);
			contig[gapDistance] = '\0';
		}
	}

	if (leftConsensusSequence != NULL && rightConsensusSequence != NULL) {
		long int extend = 10;

		gapDistance = gapDistance + 200;

		length = strlen(leftConsensusSequence);
		length1 = strlen(rightConsensusSequence);

		if (length + length1 <= gapDistance) {
			contig = (char*)malloc(sizeof(char) * (gapDistance + 1));
			strncpy(contig, leftConsensusSequence, length);
			strncpy(contig + (gapDistance - length1), rightConsensusSequence, length1);
			for (long int i = length; i < gapDistance - length1; i++) {
				contig[i] = 'N';
			}
			contig[gapDistance] = '\0';
		}
		else {
			contig = (char*)malloc(sizeof(char) * (gapDistance + extend + 1));
			long int tt = 0;
			for (long int i = 0; i < gapDistance + extend + 1; i++) {
				if (i < length) {
					contig[i] = leftConsensusSequence[i];
					tt++;
					if (tt == gapDistance) {
						for (long int p = 0; p < extend; p++) {
							contig[i + p + 1] = 'N';
						}
						break;
					}
				}
				if (i < length1) {
					contig[gapDistance + extend - i - 1] = rightConsensusSequence[length1 - i - 1];
					tt++;
					if (tt == gapDistance) {
						for (long int p = 0; p < extend; p++) {
							contig[gapDistance + extend - i - 2 - p] = 'N';
						}
						break;
					}
				}

			}
			contig[gapDistance + extend] = '\0';


		}
	}
	return contig;
}


void TransferToUpper(char* contig) {
	if (contig == NULL) {
		return;
	}
	long int contigLength = strlen(contig);
	for (long int i = 0; i < contigLength; i++) {
		if (contig[i] == 'a') {
			contig[i] = 'A';
		}
		else if (contig[i] == 't') {
			contig[i] = 'T';
		}
		else if (contig[i] == 'c') {
			contig[i] = 'C';
		}
		else if (contig[i] == 'g') {
			contig[i] = 'G';
		}
		else if (contig[i] == 'n') {
			contig[i] = 'N';
		}
	}
}

void TransferToLower(char* contig) {
	if (contig == NULL) {
		return;
	}
	long int contigLength = strlen(contig);
	for (long int i = 0; i < contigLength; i++) {
		if (contig[i] == 'A') {
			contig[i] = 'a';
		}
		else if (contig[i] == 'T') {
			contig[i] = 't';
		}
		else if (contig[i] == 'C') {
			contig[i] = 'c';
		}
		else if (contig[i] == 'G') {
			contig[i] = 'g';
		}
		else if (contig[i] == 'N') {
			contig[i] = 'n';
		}
	}
}


/*
long int GetContigIndexFromGapIndex(ScaffoldSetHead* scaffoldSetHead, long int i,long int j) {
	cout << i << "--" << j << endl;
	long int tempContigCount = 0;
	long int tempGapIndex = 0;
	long int gapLeftContigIndex = -1;
	long int gapRightContigIndex = -1;
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
		if (right == true) {//right��0
			temp = contigIndex - tempContigCount;
		}
		else {//rightΪ0
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
}*/


#endif

