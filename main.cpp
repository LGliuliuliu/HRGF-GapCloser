
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <malloc.h>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include"scaffoldSet.h"
#include "contigSet.h"
#include "aligningFromBam.h"
#include "fillGap.h"
using namespace std;


int main(int argc, char** argv) {

	if (argc == 1) {
		cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n";
		exit(0);
	}
	
	char* scaffoldSetFile = NULL;
	char* resultOutPutDirectory = NULL;
	char* hifiFile = NULL;
	long int threadCount = 1;


	int ch = 0;
	while ((ch = getopt(argc, argv, "s:r:p:t:h")) != -1) {
		switch (ch) {
		case 's': scaffoldSetFile = (char*)(optarg); break;
		case 'r': hifiFile = (char*)(optarg); break;
		case 'p': resultOutPutDirectory = (char*)optarg; break;
		case 't': threadCount = atoi(optarg); break;
		case 'h': cout << "HRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]"; exit(0); break;
		case '?': cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n"; exit(0); break;
		case ':': cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n"; exit(0); break;
		default: cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n"; exit(0); break;

		}
	}

	if (opendir(resultOutPutDirectory) == NULL) {
		mkdir(resultOutPutDirectory, 0777);
	}

	char* longScaffoldSetFile = (char*)malloc(sizeof(char) * 150);
	strcpy(longScaffoldSetFile, resultOutPutDirectory);
	strcat(longScaffoldSetFile, "/longScaffold.fa");

	char* contigSetFile = (char*)malloc(sizeof(char) * 150);
	strcpy(contigSetFile, resultOutPutDirectory);
	strcat(contigSetFile, "/contigSet.fa");

	char* aligningSam = (char*)malloc(sizeof(char) * 150);
	strcpy(aligningSam, resultOutPutDirectory);
	strcat(aligningSam, "/alingningResult.sam");
	char* aligningBam = (char*)malloc(sizeof(char) * 150);
	strcpy(aligningBam, resultOutPutDirectory);
	strcat(aligningBam, "/aligningResult.bam");

	char* file = (char*)malloc(sizeof(char) * 150);
	strcpy(file, resultOutPutDirectory);
	strcat(file, "/contigPathInHifiRead.fa");

	char* file1 = (char*)malloc(sizeof(char) * 150);
	strcpy(file1, resultOutPutDirectory);
	strcat(file1, "/optimizeContigPathInLongRead.fa");

	char* file2 = (char*)malloc(sizeof(char) * 150);
	strcpy(file2, resultOutPutDirectory);
	strcat(file2, "/supplementContigPathInLongRead.fa");

	FILE* fp;
	if ((fp = fopen(file1, "w")) == NULL) {
		printf("%s, does not exist!", file1);
		exit(0);
	}
	FILE* fp2;
	if ((fp2 = fopen(file2, "w")) == NULL) {
		printf("%s, does not exist!", file2);
		exit(0);
	}
	DelectShortContig(scaffoldSetFile, longScaffoldSetFile);//去除长度小于4000的cotnig

	//ScaffoldSetHead* scaffoldSetHead = GetScaffoldSetFromScaffoldFile(scaffoldSetFile);//获取初始scaffold以及GAP位置信息
	ScaffoldSetHead* scaffoldSetHead = GetScaffoldSetFromScaffoldFile(longScaffoldSetFile);//获取去除长度小于4000的scaffold以及GAP位置信息

	GetContigSetFromScaffoldSetHead(scaffoldSetHead);//获取contig的位置信息

	OutPutContigSetInScaffoldSetHead(scaffoldSetHead, contigSetFile);

	int maxSize = 2000;
	char* line = (char*)malloc(sizeof(char) * maxSize);

	char command[3000];

	/*sprintf(command, "minimap2 -ax asm20 %s %s > %s ", contigSetFile, hifiFile, aligningSam);
	system(command);
	sprintf(command, "samtools view -Sb %s >%s", aligningSam, aligningBam);
	system(command);*/

	long int scaffoldIndex = -1;
	long int gapIndex = -1;

	ContigSetHead* contigSetHead = GetContigSetFromContigSetFile(contigSetFile);

	GetAligningResult(scaffoldSetHead,contigSetHead, aligningBam, file);

	GetEstimatedGapLength(scaffoldSetHead, contigSetHead, file, line, maxSize);
	
	GetReadPoinstionInGap(contigSetHead, scaffoldSetHead, file, line, maxSize, fp);

	GetHiFiInGap(contigSetHead, scaffoldSetHead, hifiFile, file1, line, maxSize, resultOutPutDirectory);
	
	FillingGap(scaffoldSetHead, resultOutPutDirectory);

}