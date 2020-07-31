# SuperTAD

You may find a more detailed introduction to SuperTAD [here](https://supertad.deepomics.org)

## Install  
use git:  
```
git clone git@gitlab.deepomics.org:mengbo/SuperTAD.git SuperTAD
```
or download from source
```
wget https://supertad.deepomics.org/home/download_src SuperTAD.tar.gz
tar -xzvf SuperTAD.tar.gz
```
then
```
cd ./SuperTAD
mkdir build
cd build
cmake ..
make
```

## Dependencies 
* CMake 3.5.1+, g++, gcc

## Usage  
```
COMMANDS:
    binary  The first mode requires no user-defined parameters, run the nodes filtering by default
            ./SuperTAD <input Hi-C matrix> binary [-option values]
        OPTIONS:
            --no-filter: If given, do not filter TADs after TAD detection
    multi   The second mode requires a parameter h to determine the number of layers
            ./SuperTAD <input Hi-C matrix> multi -h <height> [-option values]
        OPTIONS:
            -h <int>: The height of coding tree, default: 2
        SHARED OPTIONS for binary and multi COMMAND:
            -K <int>: The number of leaves in the coding tree, default: nan (determined by the algorithm)
            --chrom1 <string>: chrom1 label, default: chr1
            --chrom2 <string>: chrom2 label, default: the same as chrom1
            --chrom1-start <int>: start pos on chrom1, default: 0
            --chrom2-start <int>: start pos on chrom2, default: the same as --chrom1-start
            -r/--resolution <int>: bin resolution, default: 10000
    filter  The nodes filter for optimal coding tree:
            ./SuperTAD <input Hi-C matrix> filter -i <original result> 
        OPTIONS:
            -i <string>: The list of TAD candidates
    compare The symmetric metric overlapping ratio to assess the agreement between two results
            ./SuperTAD compare <result1> <result2>
GLOBAL OPTIONS:
    -w <string>: Working directory path, default: the directory where the input Hi-C matrix is located
    -v/--verbose: Print verbose
```

## Input and Output
SuperTAD only supports dense matrix as input Hi-C matrix for now. We upload two examples of Hi-C matrix from 
[*Rao et al., Cell 2014*](https://www.cell.com/fulltext/S0092-8674(14)01497-4) as well as the results into `./data`.

The binary mode's result before filtering is stored in *.binary.original.tsv;
The binary mode's result after filtering or the filter mode's result is stored in *.binary.filter.tsv;
The multi mode's result is stored in *.multi.tsv;
All of the TAD results use the eight-column format, which records the bin indexes of detected boundaries and the genomic start and end coordinates.

An example output is shown below (resolution=1kb):
```
chr1	1	0	    1000	chr1	44	43000	44000
chr1	9	8000	9000	chr1	16	15000	16000
chr1	17	16000	17000	chr1	44	43000	44000
...
```
Each column is represented as:
1st-the chromosome of left boundary
2nd-the bin index that identified as the left boundary (start bin)
3rd-the start coordinate of start bin, in bp
4th-the end coordinate of start bin, in bp
5th-the chromosome of right boundary
6th-the bin index that identified as the right boundary (end bin)
7th-the start coordinate of end bin, in bp
8th-the end coordinate of end bin, in bp

## Examples
```
./build/SuperTAD ./data/example_sub_GM12878_chr19_KR25kb_matrix.txt binary --chrom1 chr19 -r 25000
```
This command will run binary mode (SuperTAD) on the chr19 of GM12878 at 25kb resolution and save all TADs to the example_sub_GM12878_chr19_KR25kb_matrix.txt.binary.original.tsv. 
As --no-filter is not given, the nodes filtering runs by default and save the selected TADs to the example_sub_GM12878_chr19_KR25kb_matrix.txt.binary.filter.tsv.
```
./build/SuperTAD ./data/example_sub_GM12878_chr19_KR25kb_matrix.txt multi -h 2 --chrom1 chr19 -r 25000
```
This command will run multinary mode (SuperTAD(h)) on the chr19 of GM12878 at 25kb resolution and save all TADs to the example_sub_GM12878_chr19_KR25kb_matrix.txt.multi.tsv
```
./build/SuperTAD ./data/example_sub_GM12878_chr19_KR25kb_matrix.txt filter -i ./data/example_sub_GM12878_chr19_KR25kb_matrix.txt.binary.original.tsv
```
This command will independently run the nodes filtering for the TADs in -i indicated result and save the selected TADs to *.binary.filter.tsv.
```
./build/SuperTAD compare ./data/example_sub_GM12878_chr19_KR25kb_matrix.txt.multi.tsv ./data/example_sub_IMR90_chr19_KR25kb_matrix.txt.multi.tsv
```
This command will compute the overlapping ratio between two results. 



```