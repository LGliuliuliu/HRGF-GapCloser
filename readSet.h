#ifndef READSET_H_INCLUDED 
#define READSET_H_INCLUDED 


#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

using namespace std;


typedef struct ReadSet{
	char * read;
	int readLength;	
}ReadSet;

typedef struct ReadSetHead{
	ReadSet * readSet;
	long int readCount;
}ReadSetHead;


ReadSetHead * GetReadSet(char * readSetFile,long int readCount, bool * token);

int FileIsNull(char* file);



#endif