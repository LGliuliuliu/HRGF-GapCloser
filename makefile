CC=g++ --std=c++11 -no-pie

CPPFLAGS = -g -Wall -O3

HRGF: main.o readSet.o contigSet.o scaffoldSet.o aligningFromBam.o raconConsensus.o fillGap.o
	$(CC) -o $@ $^ -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lm -ldl -lz
	
main.o: main.cpp scaffoldSet.h contigSet.h aligningFromBam.h fillGap.h
	$(CC) -c main.cpp -I $(BAMTOOLS_HOME_INCLUDE)/

readSet.o: readSet.cpp readSet.h
	$(CC) -c readSet.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

contigSet.o: contigSet.cpp contigSet.h
	$(CC) -c contigSet.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

scaffoldSet.o: scaffoldSet.cpp readSet.h contigSet.h scaffoldSet.h
	$(CC) -c scaffoldSet.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
aligningFromBam.o: aligningFromBam.cpp readSet.h contigSet.h scaffoldSet.h aligningFromBam.h
	$(CC) -c aligningFromBam.cpp -I $(BAMTOOLS_HOME_INCLUDE)/

raconConsensus.o: raconConsensus.cpp readSet.h contigSet.h scaffoldSet.h aligningFromBam.h raconConsensus.h
	$(CC) -c raconConsensus.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz
	
fillGap.o: fillGap.cpp readSet.h contigSet.h scaffoldSet.h aligningFromBam.h raconConsensus.h fillGap.h
	$(CC) -c fillGap.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ -lz

all: HRGF
	rm -f *.o

clean:
	rm -f *.o
	rm HRGF