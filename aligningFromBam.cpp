#ifndef aligningFromBam_CPP_INCLUDED 
#define aligningFromBam_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include<vector>


#include "aligningFromBam.h"

using namespace std;

int GetAligningResult(ScaffoldSetHead* scaffoldSetHead,ContigSetHead* contigSetHead, char* aligningResultBam, char* file) {//file是aligningResult.bam
	long int i = 0;
	long int j = 0;
	long int previousReadLength = 0;
	long int readLength = 0;
	long int readIndex = -1;
	long int previousReadIndex = -1;
	long int minAlignmentScore = 20;//************************************最小比对分数（参数）
	FILE* fp;
	if ((fp = fopen(file, "w")) == NULL) {
		printf("%s,does not exist!",file);
		exit(0);
	}
	string bamFileName = aligningResultBam;
	long int contigIndex = -1;
	long int previousContigIndex = -1;
	string readName;
	string previousReadName="a";
	string previousReadName1 = "a";
	BamReader bamReader;
	bamReader.Open(bamFileName);
	BamAlignment alignment;
	long int contigCount = bamReader.GetReferenceCount();//获得参考基因组的数量
	long int minContigLengthAlign = 0;//*************************************************************************比对长度（参数）
	AligningResultHead* aligningResultHead = (AligningResultHead*)malloc(sizeof(AligningResultHead));
	aligningResultHead->allocateAligningResultCount = 1000;
	aligningResultHead->aligningResultCount = 0;
	aligningResultHead->aligningShortContigResultCount = 0;
	aligningResultHead->aligningResult = (AligningResult*)malloc(sizeof(AligningResult) * aligningResultHead->allocateAligningResultCount);//构建大小为1000的结构体
	for (i = 0; i < aligningResultHead->allocateAligningResultCount; i++) {
		aligningResultHead->aligningResult[i].readStartPosition = -1;
		aligningResultHead->aligningResult[i].readEndPosition = -1;
		aligningResultHead->aligningResult[i].contigStartPosition = -1;
		aligningResultHead->aligningResult[i].contigEndPosition = -1;
		aligningResultHead->aligningResult[i].contigIndex = -1;
		aligningResultHead->aligningResult[i].readIndex = 0;
		aligningResultHead->aligningResult[i].overlapLength = -1;
		aligningResultHead->aligningResult[i].orientation = false;
	}
	i = 0;
	while (bamReader.GetNextAlignment(alignment)) {
	/*for (long int a = 0; a < 20000; a++) {
		bamReader.GetNextAlignment(alignment);*/
		readName = alignment.Name;
		//cout << readName << endl;
		contigIndex = alignment.RefID;
		/*if (contigIndex == 145) {
			cout << "-----------------------------contigIndex=" << contigIndex << endl;
		}
		if (contigIndex == 53) {
			cout << "------------------------------contigIndex=" << contigIndex << endl;
		}*/

		if (readName != previousReadName1) {//用于确定read索引
			previousReadName1 = readName;
			readIndex++;
		}
		if (!alignment.IsMapped() || alignment.MapQuality < minAlignmentScore) {//排除不正确比对
			//cout << "MapQuality=" << alignment.MapQuality << endl;
			previousReadName = readName;
			continue;
		}
		//contigIndex = alignment.RefID;
		//计算读数长度
		readLength = 0;
		const vector<CigarOp>& cigarData = alignment.CigarData;
		vector<CigarOp>::const_iterator cigarBegin = cigarData.begin();
		vector<CigarOp>::const_iterator cigarIter = cigarBegin;
		vector<CigarOp>::const_iterator cigarEnd = cigarData.end();
		for (; cigarIter != cigarEnd; ++cigarIter) {
			const CigarOp& op = (*cigarIter);
			if (op.Type != 'D') {
				readLength = readLength + op.Length;
			}
		}
		//cout << "readLength=" << readLength << endl;
		if (previousReadName != "a" && readName != previousReadName) {
			if (aligningResultHead->aligningResultCount > 0) {
				//cout << "error="<<previousReadName << endl;
				OutPutAligningResultOneLine(scaffoldSetHead,aligningResultHead, contigSetHead, fp, previousReadIndex, previousReadLength);
			}
			i = 0;
			aligningResultHead->aligningResultCount = 0;
		}
		previousReadIndex = readIndex;
		if (GetAligningResultOneLine(aligningResultHead,alignment,contigSetHead,i)!=false) {//i表示read获得比对的数量
			i++;
		}
		previousReadLength = readLength;
		previousReadName = readName;
	}
	fflush(fp);
	fclose(fp);

	return 0;


}

bool GetAligningResultOneLine(AligningResultHead* aligningResultHead, BamAlignment alignment, ContigSetHead* contigSetHead, long int index) {
	int maxAlignmentLength = 500;
	std::vector<int>clipSizes;
	std::vector<int>readPositions;
	std::vector<int>genomePositions;
	
	int readStartPosition = -1;
	int readEndPosition = -1;
	int contigStartPosition = -1;
	int contigEndPosition = -1;

	//string readName = alignment.Name;
	long int contigIndex = alignment.RefID;
	long int contigLength = contigSetHead->contigSet[contigIndex].contigLength;
	//计算读数长度
	long int readLength = 0;
	long int flag = 0;
	const vector<CigarOp>& cigarData = alignment.CigarData;
	vector<CigarOp>::const_iterator cigarBegin = cigarData.begin();
	vector<CigarOp>::const_iterator cigarIter = cigarBegin;
	vector<CigarOp>::const_iterator cigarEnd = cigarData.end();
	for (; cigarIter != cigarEnd; ++cigarIter) {
		const CigarOp& op = (*cigarIter);
		//cout << op.Type;
		
		if (op.Type != 'H' && op.Type != 'S') {
			flag = 1;
		}

		if (op.Type != 'D') {
			readLength = readLength + op.Length;
		}

		if (op.Type == 'H' || op.Type == 'S') {
			if (flag == 0) {
				readStartPosition = op.Length;
			}
			if (flag == 1) {
				readEndPosition = op.Length;
			}
		}
	}
	
	if (readStartPosition == -1) {
		readStartPosition = 0;
	}
	if (readEndPosition == -1) {
		readEndPosition = readLength - 1;
	}
	else {
		readEndPosition = readLength - readEndPosition-1;
	}
	contigStartPosition = alignment.Position;
	contigEndPosition = alignment.GetEndPosition() - 1;
	string readName = alignment.Name;


	if (readStartPosition < contigStartPosition) {
		if (readStartPosition > maxAlignmentLength) {
			//cout << "taotai1" << endl;
			return false;
		}
		aligningResultHead->aligningResult[index].readStartPosition = readStartPosition;
		aligningResultHead->aligningResult[index].contigStartPosition = contigStartPosition;
	}
	else {
		if (contigStartPosition > maxAlignmentLength) {
			//cout << "taotai2" << endl;
			return false;
		}
		aligningResultHead->aligningResult[index].readStartPosition = readStartPosition;
		aligningResultHead->aligningResult[index].contigStartPosition = contigStartPosition;
	}

	if (readLength - readEndPosition < contigLength - contigEndPosition) {
		if (readLength - readEndPosition > maxAlignmentLength) {
			//cout << "taotai3" << endl;
			return false;
		}
		aligningResultHead->aligningResult[index].readEndPosition = readEndPosition;
		aligningResultHead->aligningResult[index].contigEndPosition = contigEndPosition;
	}
	else {
		if (contigLength - contigEndPosition > maxAlignmentLength) {
			//cout << "taotai4" << endl;
			return false;
		}
		aligningResultHead->aligningResult[index].readEndPosition = readEndPosition;
		aligningResultHead->aligningResult[index].contigEndPosition = contigEndPosition ;
	}
	
	if (alignment.IsReverseStrand() != false) {//如果反向比对则得到的位置是正向比对的位置
		int temp = aligningResultHead->aligningResult[index].readStartPosition;
		aligningResultHead->aligningResult[index].readStartPosition = readLength - aligningResultHead->aligningResult[index].readEndPosition - 1;
		aligningResultHead->aligningResult[index].readEndPosition = readLength - temp - 1;
	}
	aligningResultHead->aligningResult[index].contigIndex = contigIndex;
	aligningResultHead->aligningResult[index].overlapLength = contigEndPosition - contigStartPosition + 1;
	aligningResultHead->aligningResult[index].orientation = !alignment.IsReverseStrand();

	if (aligningResultHead->aligningResult[index].overlapLength < 1000) {
		//cout << "aligningResultHead->aligningResult[index].overlapLength < 250" << endl;
		return false;
	}

	clipSizes.clear();
	readPositions.clear();
	genomePositions.clear();
	aligningResultHead->aligningResultCount++;
	return true;
}

void OutPutAligningResultOneLine(ScaffoldSetHead* scaffoldSetHead,AligningResultHead* aligningResultHead, ContigSetHead* contigSetHead, FILE* fp, long int readIndex, long int readLength) {
	int DelectNum = 0;
	
	for (long int i = 0; i < aligningResultHead->aligningResultCount - 1; i++) {//按比对开始的位置排序

		for (long int j = i + 1; j < aligningResultHead->aligningResultCount; j++) {

			if (aligningResultHead->aligningResult[i].contigIndex > aligningResultHead->aligningResult[j].contigIndex) {

				int temp = aligningResultHead->aligningResult[i].readStartPosition;
				aligningResultHead->aligningResult[i].readStartPosition = aligningResultHead->aligningResult[j].readStartPosition;
				aligningResultHead->aligningResult[j].readStartPosition = temp;

				temp = aligningResultHead->aligningResult[i].readEndPosition;
				aligningResultHead->aligningResult[i].readEndPosition = aligningResultHead->aligningResult[j].readEndPosition;
				aligningResultHead->aligningResult[j].readEndPosition = temp;

				temp = aligningResultHead->aligningResult[i].contigStartPosition;
				aligningResultHead->aligningResult[i].contigStartPosition = aligningResultHead->aligningResult[j].contigStartPosition;
				aligningResultHead->aligningResult[j].contigStartPosition = temp;

				temp = aligningResultHead->aligningResult[i].contigEndPosition;
				aligningResultHead->aligningResult[i].contigEndPosition = aligningResultHead->aligningResult[j].contigEndPosition;
				aligningResultHead->aligningResult[j].contigEndPosition = temp;

				temp = aligningResultHead->aligningResult[i].contigIndex;
				aligningResultHead->aligningResult[i].contigIndex = aligningResultHead->aligningResult[j].contigIndex;
				aligningResultHead->aligningResult[j].contigIndex = temp;

				temp = aligningResultHead->aligningResult[i].overlapLength;
				aligningResultHead->aligningResult[i].overlapLength = aligningResultHead->aligningResult[j].overlapLength;
				aligningResultHead->aligningResult[j].overlapLength = temp;

				temp = aligningResultHead->aligningResult[i].orientation;
				aligningResultHead->aligningResult[i].orientation = aligningResultHead->aligningResult[j].orientation;
				aligningResultHead->aligningResult[j].orientation = temp;
			}
			if (aligningResultHead->aligningResult[i].contigIndex == aligningResultHead->aligningResult[j].contigIndex && aligningResultHead->aligningResult[i].contigStartPosition > aligningResultHead->aligningResult[j].contigStartPosition) {
				int temp = aligningResultHead->aligningResult[i].readStartPosition;
				aligningResultHead->aligningResult[i].readStartPosition = aligningResultHead->aligningResult[j].readStartPosition;
				aligningResultHead->aligningResult[j].readStartPosition = temp;

				temp = aligningResultHead->aligningResult[i].readEndPosition;
				aligningResultHead->aligningResult[i].readEndPosition = aligningResultHead->aligningResult[j].readEndPosition;
				aligningResultHead->aligningResult[j].readEndPosition = temp;

				temp = aligningResultHead->aligningResult[i].contigStartPosition;
				aligningResultHead->aligningResult[i].contigStartPosition = aligningResultHead->aligningResult[j].contigStartPosition;
				aligningResultHead->aligningResult[j].contigStartPosition = temp;

				temp = aligningResultHead->aligningResult[i].contigEndPosition;
				aligningResultHead->aligningResult[i].contigEndPosition = aligningResultHead->aligningResult[j].contigEndPosition;
				aligningResultHead->aligningResult[j].contigEndPosition = temp;

				temp = aligningResultHead->aligningResult[i].contigIndex;
				aligningResultHead->aligningResult[i].contigIndex = aligningResultHead->aligningResult[j].contigIndex;
				aligningResultHead->aligningResult[j].contigIndex = temp;

				temp = aligningResultHead->aligningResult[i].overlapLength;
				aligningResultHead->aligningResult[i].overlapLength = aligningResultHead->aligningResult[j].overlapLength;
				aligningResultHead->aligningResult[j].overlapLength = temp;

				temp = aligningResultHead->aligningResult[i].orientation;
				aligningResultHead->aligningResult[i].orientation = aligningResultHead->aligningResult[j].orientation;
				aligningResultHead->aligningResult[j].orientation = temp;
			}
		}
	}

	
	
	//for (long int i = 0; i < aligningResultHead->aligningResultCount - 1; i++) {//同一个contig比对到同一个读数的两个位置则去掉该比对

	//	for (long int j = i + 1; j < aligningResultHead->aligningResultCount; j++) {

	//		if (aligningResultHead->aligningResult[i].contigIndex == aligningResultHead->aligningResult[j].contigIndex && aligningResultHead->aligningResult[i].contigIndex != -1) {

	//			//cout << aligningResultHead->aligningResult[i].contigIndex << "---0---" << aligningResultHead->aligningResult[i].orientation << endl;
	//			//cout << i << "    " << j << endl;
	//			/*cout << "reaptaligncontig=" << aligningResultHead->aligningResult[i].contigIndex << endl;
	//			cout << "004------contigStart=" << aligningResultHead->aligningResult[i].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[i].contigEndPosition << endl;
	//			cout << "004------contigStart=" << aligningResultHead->aligningResult[j].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[j].contigEndPosition << endl;*/

	//			if ((aligningResultHead->aligningResult[i].contigStartPosition < aligningResultHead->aligningResult[j].contigStartPosition && aligningResultHead->aligningResult[i].contigEndPosition > aligningResultHead->aligningResult[j].contigEndPosition)
	//				|| (aligningResultHead->aligningResult[i].contigStartPosition > aligningResultHead->aligningResult[j].contigStartPosition && aligningResultHead->aligningResult[i].contigEndPosition < aligningResultHead->aligningResult[j].contigEndPosition)) {
	//				/*cout << "001------contigStart=" << aligningResultHead->aligningResult[i].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[i].contigEndPosition << endl;
	//				cout << "001------contigStart=" << aligningResultHead->aligningResult[j].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[j].contigEndPosition << endl;*/

	//				if (aligningResultHead->aligningResult[i].overlapLength > aligningResultHead->aligningResult[j].overlapLength) {

	//					aligningResultHead->aligningResult[j].contigIndex = -1;
	//					DelectNum = DelectNum + 1;
	//					//cout << "   1" << endl;
	//					
	//				}
	//				else {
	//					aligningResultHead->aligningResult[i].contigIndex = -1;
	//					DelectNum = DelectNum + 1;
	//					//cout << "   2" << endl;
	//					
	//				}

	//			}
	//			if ((aligningResultHead->aligningResult[i].contigStartPosition == aligningResultHead->aligningResult[j].contigStartPosition && aligningResultHead->aligningResult[i].contigEndPosition != aligningResultHead->aligningResult[j].contigEndPosition)
	//				|| (aligningResultHead->aligningResult[i].contigStartPosition != aligningResultHead->aligningResult[j].contigStartPosition && aligningResultHead->aligningResult[i].contigEndPosition == aligningResultHead->aligningResult[j].contigEndPosition)) {
	//				/*cout << "002------contigStart=" << aligningResultHead->aligningResult[i].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[i].contigEndPosition << endl;
	//				cout << "002------contigStart=" << aligningResultHead->aligningResult[j].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[j].contigEndPosition << endl;*/

	//				if (aligningResultHead->aligningResult[i].overlapLength > aligningResultHead->aligningResult[j].overlapLength) {

	//					aligningResultHead->aligningResult[j].contigIndex = -1;
	//					DelectNum = DelectNum + 1;
	//					//cout << "   3" << endl;
	//					
	//				}
	//				else {

	//					aligningResultHead->aligningResult[i].contigIndex = -1;
	//					DelectNum = DelectNum + 1;
	//					//cout << "   4" << endl;
	//					
	//				}
	//			}
	//		}
	//	}
	//}



	//for (long int i = 0; i < aligningResultHead->aligningResultCount - 1; i++) {//同一个contig比对到同一个读数的两个位置则去掉该比对

	//	for (long int j = i + 1; j < aligningResultHead->aligningResultCount; j++) {

	//		if (aligningResultHead->aligningResult[i].contigIndex == aligningResultHead->aligningResult[j].contigIndex && aligningResultHead->aligningResult[i].contigIndex != -1) {

	//			//cout << aligningResultHead->aligningResult[i].contigIndex << "---2---" << aligningResultHead->aligningResult[i].orientation << endl;
	//			//cout << i << "    " << j << endl;
	//			//cout << "reaptaligncontig=" << aligningResultHead->aligningResult[i].contigIndex << endl;
	//			if (aligningResultHead->aligningResult[i].contigStartPosition == aligningResultHead->aligningResult[j].contigStartPosition && aligningResultHead->aligningResult[i].contigEndPosition == aligningResultHead->aligningResult[j].contigEndPosition) {
	//				/*cout << "004------contigStart=" << aligningResultHead->aligningResult[i].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[i].contigEndPosition << endl;
	//				cout << "004------contigStart=" << aligningResultHead->aligningResult[j].contigStartPosition << "     contigEnd=" << aligningResultHead->aligningResult[j].contigEndPosition << endl;*/
	//				
	//				//通过gap长度排除相同contig
	//				for (long int k = 0; k < aligningResultHead->aligningResultCount; k++) {

	//					long int scaffoldIndex = -1;
	//					long int gapIndex = -1;

	//					if ((aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[i].contigIndex - 1 || aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[i].contigIndex + 1)
	//						&& aligningResultHead->aligningResult[k].orientation == aligningResultHead->aligningResult[i].orientation) {

	//						//cout << aligningResultHead->aligningResult[k].contigIndex << "-----" << aligningResultHead->aligningResult[i].contigIndex << endl;

	//						GetGapIndexFromContigIndex(scaffoldSetHead, aligningResultHead->aligningResult[k].contigIndex, aligningResultHead->aligningResult[i].contigIndex, scaffoldIndex, gapIndex);

	//						long int gapLength = 0;
	//						long int realGapLength01 = 0;
	//						long int realGapLength02 = 0;
	//						long int diff01 = 0;
	//						long int diff02 = 0;
	//						//cout << scaffoldIndex << "**************" << gapIndex << endl;

	//						if (scaffoldIndex != -1 && gapIndex != -1) {
	//							gapLength = scaffoldSetHead->scaffoldSet[scaffoldIndex].gapSet[gapIndex].gapLength;
	//							/*cout << "gapLength=" << gapLength << endl;
	//							cout << "k:         " << aligningResultHead->aligningResult[k].readStartPosition << "------------------" << aligningResultHead->aligningResult[k].readEndPosition << endl;
	//							cout << "j:          " << aligningResultHead->aligningResult[j].readStartPosition << "------------------" << aligningResultHead->aligningResult[j].readEndPosition << endl;
	//							cout << "i:          " << aligningResultHead->aligningResult[i].readStartPosition << "------------------" << aligningResultHead->aligningResult[i].readEndPosition << endl;*/


	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[i].contigIndex - 1 && aligningResultHead->aligningResult[i].orientation==1 
	//								&& aligningResultHead->aligningResult[k].readStartPosition < aligningResultHead->aligningResult[i].readStartPosition) {
	//								realGapLength01 = aligningResultHead->aligningResult[i].readStartPosition - aligningResultHead->aligningResult[k].readEndPosition;
	//							}

	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[i].contigIndex - 1 && aligningResultHead->aligningResult[i].orientation == 0 
	//								&& aligningResultHead->aligningResult[i].readStartPosition < aligningResultHead->aligningResult[k].readStartPosition) {
	//								realGapLength01 = aligningResultHead->aligningResult[k].readStartPosition - aligningResultHead->aligningResult[i].readEndPosition;
	//							}

	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[i].contigIndex + 1 && aligningResultHead->aligningResult[i].orientation == 1
	//								&& aligningResultHead->aligningResult[k].readStartPosition > aligningResultHead->aligningResult[i].readStartPosition) {
	//								realGapLength01 = aligningResultHead->aligningResult[k].readStartPosition - aligningResultHead->aligningResult[i].readEndPosition;
	//							}

	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[i].contigIndex + 1 && aligningResultHead->aligningResult[i].orientation == 0
	//								&& aligningResultHead->aligningResult[i].readStartPosition > aligningResultHead->aligningResult[k].readStartPosition) {
	//								realGapLength01 = aligningResultHead->aligningResult[i].readStartPosition - aligningResultHead->aligningResult[k].readEndPosition;
	//							}

	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[j].contigIndex - 1 && aligningResultHead->aligningResult[j].orientation == 1
	//								&& aligningResultHead->aligningResult[k].readStartPosition < aligningResultHead->aligningResult[j].readStartPosition) {
	//								realGapLength02 = aligningResultHead->aligningResult[j].readStartPosition - aligningResultHead->aligningResult[k].readEndPosition;
	//							}

	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[j].contigIndex - 1 && aligningResultHead->aligningResult[j].orientation == 0
	//								&& aligningResultHead->aligningResult[k].readStartPosition > aligningResultHead->aligningResult[j].readStartPosition) {
	//								realGapLength02 = aligningResultHead->aligningResult[k].readStartPosition - aligningResultHead->aligningResult[j].readEndPosition;
	//							}

	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[j].contigIndex + 1 && aligningResultHead->aligningResult[j].orientation == 1
	//								&& aligningResultHead->aligningResult[k].readStartPosition > aligningResultHead->aligningResult[j].readStartPosition) {
	//								realGapLength02 = aligningResultHead->aligningResult[k].readStartPosition - aligningResultHead->aligningResult[j].readEndPosition;
	//							}

	//							if (aligningResultHead->aligningResult[k].contigIndex == aligningResultHead->aligningResult[j].contigIndex + 1 && aligningResultHead->aligningResult[j].orientation == 0
	//								&& aligningResultHead->aligningResult[j].readStartPosition > aligningResultHead->aligningResult[k].readStartPosition) {
	//								realGapLength02 = aligningResultHead->aligningResult[j].readStartPosition - aligningResultHead->aligningResult[k].readEndPosition;
	//							}
	//							//cout << realGapLength01 << "===========" << realGapLength02 << endl;
	//							diff01 = abs(realGapLength01 - gapLength);
	//							diff02 = abs(realGapLength02 - gapLength);
	//							//cout << diff01 << "==========" << diff02 << endl;

	//							if (diff01 > diff02) {
	//								aligningResultHead->aligningResult[i].contigIndex = -1;
	//								DelectNum = DelectNum + 1;
	//								//cout << "   5" << endl;
	//							}if (diff01 < diff02) {
	//								aligningResultHead->aligningResult[j].contigIndex = -1;
	//								DelectNum = DelectNum + 1;
	//								//cout << "   6" << endl;
	//							}
	//							if (diff01 == diff02) {
	//								aligningResultHead->aligningResult[i].contigIndex = -1;
	//								aligningResultHead->aligningResult[j].contigIndex = -1;
	//								DelectNum = DelectNum + 2;
	//								//cout << "   7" << endl;
	//							}

	//						}
	//					}
	//				}
	//				//如果排除不成功则去掉同一个contig的不同比对
	//				if (aligningResultHead->aligningResult[i].contigIndex == aligningResultHead->aligningResult[j].contigIndex) {
	//					aligningResultHead->aligningResult[i].contigIndex = -1;
	//					aligningResultHead->aligningResult[j].contigIndex = -1;
	//					DelectNum = DelectNum + 2;
	//					//cout << "   8" << endl;
	//				}
	//			}

	//			//合并相同两个contig的两个比对
	//			if (aligningResultHead->aligningResult[i].orientation == aligningResultHead->aligningResult[j].orientation && aligningResultHead->aligningResult[i].orientation == 1
	//				&& aligningResultHead->aligningResult[i].contigStartPosition < aligningResultHead->aligningResult[j].contigStartPosition && aligningResultHead->aligningResult[i].contigEndPosition < aligningResultHead->aligningResult[j].contigEndPosition
	//				&& aligningResultHead->aligningResult[i].readStartPosition < aligningResultHead->aligningResult[j].readStartPosition && aligningResultHead->aligningResult[i].readEndPosition < aligningResultHead->aligningResult[j].readEndPosition) {

	//				/*cout << aligningResultHead->aligningResult[i].contigIndex << "----3--" << aligningResultHead->aligningResult[i].orientation << "---------" << aligningResultHead->aligningResult[j].orientation << endl;
	//				cout << aligningResultHead->aligningResult[i].contigStartPosition << "--------------" << aligningResultHead->aligningResult[i].contigEndPosition << "-----------" << aligningResultHead->aligningResult[i].readStartPosition << "-----------" << aligningResultHead->aligningResult[i].readEndPosition << endl;
	//				cout << aligningResultHead->aligningResult[j].contigStartPosition << "--------------" << aligningResultHead->aligningResult[j].contigEndPosition << "-----------" << aligningResultHead->aligningResult[j].readStartPosition << "-----------" << aligningResultHead->aligningResult[j].readEndPosition << endl;*/


	//				aligningResultHead->aligningResult[i].contigEndPosition = aligningResultHead->aligningResult[j].contigEndPosition;
	//				aligningResultHead->aligningResult[i].readEndPosition = aligningResultHead->aligningResult[j].readEndPosition;
	//				aligningResultHead->aligningResult[i].overlapLength = aligningResultHead->aligningResult[i].contigEndPosition - aligningResultHead->aligningResult[i].contigStartPosition;
	//				aligningResultHead->aligningResult[j].contigIndex = -1;
	//				DelectNum = DelectNum + 1;
	//				//cout << "   9" << endl;
	//				/*cout << "newContigStart=" << aligningResultHead->aligningResult[i].contigStartPosition << "     newcontigEnd=" << aligningResultHead->aligningResult[i].contigEndPosition << endl;
	//				cout << "newReadStart=" << aligningResultHead->aligningResult[i].readStartPosition << "     newReadEnd=" << aligningResultHead->aligningResult[i].readEndPosition << endl;
	//				cout << "newOverlap=" << aligningResultHead->aligningResult[i].overlapLength << endl;*/
	//			}
	//			//合并相同两个contig的两个比对
	//			if (aligningResultHead->aligningResult[i].orientation == aligningResultHead->aligningResult[j].orientation && aligningResultHead->aligningResult[i].orientation == 0
	//				&& aligningResultHead->aligningResult[i].contigStartPosition < aligningResultHead->aligningResult[j].contigStartPosition && aligningResultHead->aligningResult[i].contigEndPosition < aligningResultHead->aligningResult[j].contigEndPosition
	//				&& aligningResultHead->aligningResult[i].readStartPosition > aligningResultHead->aligningResult[j].readStartPosition && aligningResultHead->aligningResult[i].readEndPosition > aligningResultHead->aligningResult[j].readEndPosition) {

	//				/*cout << aligningResultHead->aligningResult[i].contigIndex << "----3--" << aligningResultHead->aligningResult[i].orientation << "---------" << aligningResultHead->aligningResult[j].orientation << endl;
	//				cout << aligningResultHead->aligningResult[i].contigStartPosition << "--------------" << aligningResultHead->aligningResult[i].contigEndPosition << "-----------" << aligningResultHead->aligningResult[i].readStartPosition << "-----------" << aligningResultHead->aligningResult[i].readEndPosition << endl;
	//				cout << aligningResultHead->aligningResult[j].contigStartPosition << "--------------" << aligningResultHead->aligningResult[j].contigEndPosition << "-----------" << aligningResultHead->aligningResult[j].readStartPosition << "-----------" << aligningResultHead->aligningResult[j].readEndPosition << endl;*/

	//				aligningResultHead->aligningResult[i].contigEndPosition = aligningResultHead->aligningResult[j].contigEndPosition;
	//				aligningResultHead->aligningResult[i].readStartPosition = aligningResultHead->aligningResult[j].readStartPosition;
	//				aligningResultHead->aligningResult[i].overlapLength = aligningResultHead->aligningResult[i].contigEndPosition - aligningResultHead->aligningResult[i].contigStartPosition;
	//				aligningResultHead->aligningResult[j].contigIndex = -1;
	//				DelectNum = DelectNum + 1;
	//				//cout << "   10" << endl;
	//				/*cout << "newContigStart1=" << aligningResultHead->aligningResult[i].contigStartPosition << "     newcontigEnd1=" << aligningResultHead->aligningResult[i].contigEndPosition << endl;
	//				cout << "newReadStart1=" << aligningResultHead->aligningResult[i].readStartPosition << "     newReadEnd1=" << aligningResultHead->aligningResult[i].readEndPosition << endl;
	//				cout << "newOverlap1=" << aligningResultHead->aligningResult[i].overlapLength << endl;*/
	//			}


	//			
	//		}
	//	}

	//}

	//for (long int i = 0; i < aligningResultHead->aligningResultCount - 1; i++) {//同一个contig比对到同一个读数的两个位置则去掉该比对

	//	for (long int j = i + 1; j < aligningResultHead->aligningResultCount; j++) {
	//		if (aligningResultHead->aligningResult[i].contigIndex == aligningResultHead->aligningResult[j].contigIndex && aligningResultHead->aligningResult[i].contigIndex != -1) {
	//			cout << "readIndex=" << readIndex << endl;
	//			cout << "i=" << i << "----j=" << j << endl;
	//			cout << "contigLength=" << contigSetHead->contigSet[aligningResultHead->aligningResult[i].contigIndex].contigLength << "  readLength=" << readLength << endl;
	//			cout << "contigIndex=" << aligningResultHead->aligningResult[i].contigIndex << "  orine:" << aligningResultHead->aligningResult[i].orientation << "---------" << aligningResultHead->aligningResult[j].orientation << endl;
	//			cout << "i:" << aligningResultHead->aligningResult[i].contigStartPosition << "------" << aligningResultHead->aligningResult[i].contigEndPosition
	//				<< "---------" << aligningResultHead->aligningResult[i].readStartPosition << "---------" << aligningResultHead->aligningResult[i].readEndPosition << "---------" << aligningResultHead->aligningResult[i].overlapLength << endl;

	//			cout << "j:" << aligningResultHead->aligningResult[j].contigStartPosition << "------" << aligningResultHead->aligningResult[j].contigEndPosition
	//				<< "---------" << aligningResultHead->aligningResult[j].readStartPosition << "---------" << aligningResultHead->aligningResult[j].readEndPosition << "---------" << aligningResultHead->aligningResult[j].overlapLength << endl;
	//			cout << endl;
	//		}
	//	}
	//}


	
	for (long int i = 0; i < aligningResultHead->aligningResultCount - 1; i++) {//去除两个contig同时比对到读数头或尾

		for (long int j = i + 1; j < aligningResultHead->aligningResultCount; j++) {

			if (aligningResultHead->aligningResult[i].readStartPosition == aligningResultHead->aligningResult[j].readStartPosition && aligningResultHead->aligningResult[i].readStartPosition == 0 && aligningResultHead->aligningResult[i].contigIndex != -1 && aligningResultHead->aligningResult[j].contigIndex != -1) {

				if (aligningResultHead->aligningResult[i].overlapLength > aligningResultHead->aligningResult[j].overlapLength) {
					aligningResultHead->aligningResult[j].contigIndex = -1;
					DelectNum = DelectNum + 1;

				}
				if (aligningResultHead->aligningResult[i].overlapLength < aligningResultHead->aligningResult[j].overlapLength) {
					aligningResultHead->aligningResult[i].contigIndex = -1;
					DelectNum = DelectNum + 1;

				}
				if (aligningResultHead->aligningResult[i].overlapLength == aligningResultHead->aligningResult[j].overlapLength) {
					aligningResultHead->aligningResult[i].contigIndex = -1;
					aligningResultHead->aligningResult[j].contigIndex = -1;
					DelectNum = DelectNum + 2;

				}
			}
			if (aligningResultHead->aligningResult[i].readEndPosition == aligningResultHead->aligningResult[j].readEndPosition && aligningResultHead->aligningResult[i].readEndPosition == readLength - 1 && aligningResultHead->aligningResult[i].contigIndex != -1 && aligningResultHead->aligningResult[j].contigIndex != -1) {

				if (aligningResultHead->aligningResult[i].overlapLength > aligningResultHead->aligningResult[j].overlapLength) {
					aligningResultHead->aligningResult[j].contigIndex = -1;
					DelectNum = DelectNum + 1;

				}
				if (aligningResultHead->aligningResult[i].overlapLength < aligningResultHead->aligningResult[j].overlapLength) {
					aligningResultHead->aligningResult[i].contigIndex = -1;
					DelectNum = DelectNum + 1;

				}
				if (aligningResultHead->aligningResult[i].overlapLength == aligningResultHead->aligningResult[j].overlapLength) {
					aligningResultHead->aligningResult[i].contigIndex = -1;
					aligningResultHead->aligningResult[j].contigIndex = -1;
					DelectNum = DelectNum + 2;

				}
			}
		}

	}


	
	if (aligningResultHead->aligningResultCount - DelectNum > 0) {
		fprintf(fp, "%d,%ld,%ld", aligningResultHead->aligningResultCount - DelectNum, readIndex, readLength);
	}

	for (long int i = 0; i < aligningResultHead->aligningResultCount; i++) {

		if (aligningResultHead->aligningResult[i].contigIndex != -1) {
			fprintf(fp, ",%d,%d,%d,%d,%d,%d,%d", aligningResultHead->aligningResult[i].contigIndex, aligningResultHead->aligningResult[i].readStartPosition,
				aligningResultHead->aligningResult[i].readEndPosition, aligningResultHead->aligningResult[i].contigStartPosition,
				aligningResultHead->aligningResult[i].contigEndPosition, aligningResultHead->aligningResult[i].overlapLength,
				aligningResultHead->aligningResult[i].orientation);
		}

	}
	fprintf(fp, ",\n");
	
	
}


















#endif
