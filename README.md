# HRGF
HRGF is a gap filling method based on HiFi read and read clustering
=========


gap-filling: HRGF
=================

1) Introduction
```
    HRGF is an gap-filling tool which aims to filling gaps of scaffolds. 
    The input data of HRGF is the HiFi reads (fasta format) and the scaffolds (fasta format). 
```
2) Before installing and running
```
    Please install minimap2 from https://github.com/lh3/minimap2.
	Please install Samtools from https://sourceforge.net/projects/samtools/files/samtools/.
	Please build and install Bamtools from https://github.com/pezmaster31/bamtools.
	Please install blast from https://github.com/julianshapiro/blast.
```
3) Installing.
```
    HRGF should run on Linux operating sysetm with gcc. We test HRGF using gcc9.4.0 on Ubuntu.
    Create a main directory. Copy all source code to this directory.
	cd HRGF
	export BAMTOOLS_HOME_INCLUDE=/path_bamtools_include_api_shared/
	export BAMTOOLS_HOME_LIB=/path_bamtools_lib_libbamtools.a/
	make all
```
4) Running.
```
	HRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]
    
    -s <scaffolds.fa>: 
	    The file includes scaffolds produced by one assembler.
	-r <readType>: 
	    The type of HiFi reads.
	-p <outdir>: 
	   The output path of the file.
	-t <threadNumber>: 
	   The number of threads used by the tool.
```
5) Output.
```
    The output file "gapfilling-scaffoldset.fa" is the gap-filling result. 
```


