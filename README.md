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
    -w <string>: Working directory path, default: the directory where the input file is located
    -v/--verbose: Print verbose
```
The binary mode's result before filtering is stored in *.binary.original.tsv
The binary mode's result after filtering or the filter mode's result are stored in *.binary.filter.tsv
The multi mode's result is stored in *.multi.tsv

## Sample data
We only supports dense matrix as input for now.
A simulated contact matrix and a matrix derived from a partition of chromosome subtracted from 
[*Rao et al, Cell 2014*](https://www.cell.com/fulltext/S0092-8674(14)01497-4)
are included in `./data`.

For example,
By running SuperTAD on `./data/sub_chr6_200_KR100kb_matrix.txt`, 
you first need to compile the program following the first section.
Then pass the path of the input along with other parameters to the executable.

The basic command can be like this:
```
./build/SuperTAD ./data/sub_chr6_200_KR100kb_matrix.txt -b -r 1000
```

If you are not providing working directory, SuperTAD will take the directory where input lies in.

## Interpret output
SuperTAD provides a uniform eight-column format when analyzing and comparing all the methods. 
The format kept the bin coordinates of identified boundaries and the bin resolution from contact map input and each 
column is separated by tab. 

A sample output looks like below:
```
chr1	1	0	1000	chr1	44	43000	44000
chr1	1	0	1000	chr1	16	15000	16000
chr1	1	0	1000	chr1	8	7000	8000
chr1	9	8000	9000	chr1	16	15000	16000
chr1	17	16000	17000	chr1	44	43000	44000
chr1	17	16000	17000	chr1	35	34000	35000
chr1	17	16000	17000	chr1	25	24000	25000
chr1	26	25000	26000	chr1	35	34000	35000
chr1	36	35000	36000	chr1	44	43000	44000
chr1	45	44000	45000	chr1	100	99000	100000
chr1	45	44000	45000	chr1	60	59000	60000
chr1	45	44000	45000	chr1	52	51000	52000
chr1	53	52000	53000	chr1	60	59000	60000
chr1	61	60000	61000	chr1	100	99000	100000
chr1	61	60000	61000	chr1	79	78000	79000
chr1	61	60000	61000	chr1	66	65000	66000
chr1	67	66000	67000	chr1	79	78000	79000
chr1	67	66000	67000	chr1	71	70000	71000
chr1	72	71000	72000	chr1	79	78000	79000
chr1	80	79000	80000	chr1	100	99000	100000
chr1	80	79000	80000	chr1	86	85000	86000
chr1	87	86000	87000	chr1	100	99000	100000
chr1	87	86000	87000	chr1	95	94000	95000
chr1	96	95000	96000	chr1	100	99000	100000
```
The output file will has either `.binary.original.tsv` or `.binary.filter.tsv` as post-fix for binary mode, 
and `.multi.original.tsv` for multi-nary mode.

#
You may find more detailed introduction to SuperTAD [here](https://supertad.deepomics.org)