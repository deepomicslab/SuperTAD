# SuperTAD

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
mkdir build; cd build
cmake ..
make
```

## Dependencies 
* CMake 3.5.1+, g++, gcc
* [Eigen3](https://eigen.tuxfamily.org/dox/) (included in the source)

## Basic Usage  
```
usage: ./SuperTAD <input hic matrix> [options (values)]
OPTIONS:
    -w <string>: working directory; If not given, use current input directory)\n
    -b: binary tree version\n
    -m: multiple tree version\n
    -H: multiple version (h1 fast mode)\n
    -k <int>: number of leaves in candidate coding tree (default NAN)\n
    -h <int>: hierarchy number (default 2)\n
    --no-filter: do not filter TADs\n
    --no-bold: disable bold mode\n
    --bedpe: write output in BEDPE format\n
    --chrom1 <string>: chrom1 label\n
    --chrom2 <string>: chrom2 label (if only chrom1 is given, assume chrom1 and 2 are identical)\n
    --chrom1-start <int>: start pos on chrom1\n
    --chrom2-start <int>: start pos on chrom2\n
    -r/--resolution <int>: resolution\n
    -v/--verbose: print verbose\n
```

## Sample data
sample input data can be found in data (data/simulate_matrix.txt).

## Interpret output
sample output in plain text format:
```
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 
12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 
12 13 14 15 16 17 18 19 20 21 22 23 24 
36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 
36 37 38 39 40 41 42 43 44 45 46 47 48 
49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 
49 50 51 52 53 54 55 56 57 58 59 60 
```
each line represents one TAD of which boundary starts with the bin of the first index and ends with the bin of the last index in that line.

sample output in BEDPE foramt:
```
#chrom1	start1	end1	chrom2	start2	end2	name
#SuperTAD v1.1
1	1	35	1	1	35	node1
1	12	35	1	12	35	node2
1	12	24	1	12	24	node3
1	36	70	1	36	70	node4
1	36	48	1	36	48	node5
1	49	70	1	49	70	node6
1	49	60	1	49	60	node7
```

#
You may find detailed introduction to SuperTAD [here](https://supertad.deepomics.org)