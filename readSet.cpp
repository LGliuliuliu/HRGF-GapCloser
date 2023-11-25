#ifndef READSET_CPP_INCLUDED 
#define READSET_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "readSet.h"

using namespace std;

ReadSetHead * GetReadSet(char * readSetFile, long int readCount, bool * token){
    ReadSetHead* readSetHead = (ReadSetHead*)malloc(sizeof(ReadSetHead));
    readSetHead->readSet = NULL;
    readSetHead->readCount = readCount;
    readSetHead->readSet = (ReadSet*)malloc(sizeof(ReadSet) * readSetHead->readCount);
    for (long int i = 0; i < readSetHead->readCount; i++) {
        readSetHead->readSet[i].read = NULL;
        readSetHead->readSet[i].readLength = 0;
    }
    long int maxSize = 90000;
    char* read = NULL;
    if (NULL == (read = (char*)malloc(sizeof(char) * maxSize))) {
        perror("malloc error!");
        exit(1);
    }
    FILE* fp;
    if ((fp = fopen(readSetFile, "r")) == NULL) {
        printf("%s, does not exist!", readSetFile);
        exit(0);
    }
    long int allocateLength = 0;
    long int readIndex = -1;
    while ((fgets(read, maxSize, fp)) != NULL) {
        if (read[0] == '>') {
            readIndex++;
            continue;
        }
        if (token[readIndex] != true) {//如果不是需要的读数
            continue;
        }
        long int extendLength = strlen(read);
        if (read[extendLength - 1] == '\n') {
            extendLength--;
        }
        long int readLength = 0;

        char* tempRead = NULL;
        if (readSetHead->readSet[readIndex].read != NULL) {
            if (readSetHead->readSet[readIndex].readLength + extendLength >= allocateLength) {
                readLength = readSetHead->readSet[readIndex].readLength;
                readSetHead->readSet[readIndex].read = (char*)realloc(readSetHead->readSet[readIndex].read, allocateLength + maxSize + 1);

                allocateLength = allocateLength + maxSize + 1;

                strncpy(readSetHead->readSet[readIndex].read + readLength, read, extendLength);
                readSetHead->readSet[readIndex].read[readLength + extendLength] = '\0';
                readSetHead->readSet[readIndex].readLength = readLength + extendLength;

            }
            else {
                strncpy(readSetHead->readSet[readIndex].read + readSetHead->readSet[readIndex].readLength, read, extendLength);
                readSetHead->readSet[readIndex].read[readSetHead->readSet[readIndex].readLength + extendLength] = '\0';
                readSetHead->readSet[readIndex].readLength = readSetHead->readSet[readIndex].readLength + extendLength;
            }

        }
        else {
            readSetHead->readSet[readIndex].read = (char*)malloc(sizeof(char) * (maxSize + 1));
            strncpy(readSetHead->readSet[readIndex].read, read, extendLength);
            readSetHead->readSet[readIndex].read[extendLength] = '\0';
            readSetHead->readSet[readIndex].readLength = extendLength;
            allocateLength = maxSize + 1;
        }
    }
    fflush(fp);
    fclose(fp);
    
    //for (int k = 0; k < 200; k++) {
    //    if (token[k] == true) {
    //        cout << k << endl;
    //        /*cout << readSetHead->readSet[k].read << endl;*/
    //    }
    //    
    //}
    return readSetHead;

}

int FileIsNull(char* file) {
    FILE* fp = fopen(file, "r");
    char ch = fgetc(fp);
    fclose(fp);
    if (ch == EOF) {
        return 1;
    }
    else {
        return 0;
    }


}





#endif