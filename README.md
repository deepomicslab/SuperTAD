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

## Usage  
```
usage: SuperTAD <input hic matrix> [options (values)]
OPTIONS:
    -w <string>: working directory; If not given, use current input directory);
    -b: binary tree version;
    -m: multiple tree version;
    -H: multiple version (h1 fast mode);
    -k <int>: number of leaves in candidate coding tree (default NAN);
    -h <int>: hierarchy number (default 2);
    --no-filter: do not filter TADs;
    --no-fast: disable fast mode for binary mode;
    --chrom1 <string>: chrom1 label;
    --chrom2 <string>: chrom2 label (if only chrom1 is given, assume chrom1 and 2 are identical);
    --chrom1-start <int>: start pos on chrom1;
    --chrom2-start <int>: start pos on chrom2;
    -r/--resolution <int>: resolution;
    -v/--verbose: print verbose;
```

## Sample data
We only supports dense matrix as input for now.
A simulated contact matrix and a matrix derived from a partition of chromosome subtracted from [*Rao et al, Cell 2014*](https://www.cell.com/fulltext/S0092-8674(14)01497-4)
are included in ./data .

For example,
By running SuperTAD on ./data/sub_chr6_200_KR100kb_matrix.txt, 
you first need to compile the program following the first section.
Then pass the path of the input along with other parameters to the executable.

For simplicity, the basic command looks like this:

./build/SuperTAD ./data/sub_chr6_200_KR100kb_matrix.txt -b -r 1000

If you are not providing working directory, SuperTAD will take the directory where input lies in.

## Interpret output
Sample output (with sub_chr6_200_KR100kb_matrix.txt as input, with parameters '-r 1000'):
```
1	0	0	1000	43	43000	44000
1	0	0	1000	15	15000	16000
1	0	0	1000	7	7000	8000
1	8	8000	9000	15	15000	16000
1	16	16000	17000	43	43000	44000
1	16	16000	17000	34	34000	35000
1	16	16000	17000	24	24000	25000
1	25	25000	26000	34	34000	35000
1	35	35000	36000	43	43000	44000
1	44	44000	45000	99	99000	100000
1	44	44000	45000	59	59000	60000
1	44	44000	45000	51	51000	52000
1	52	52000	53000	59	59000	60000
1	60	60000	61000	99	99000	100000
1	60	60000	61000	78	78000	79000
1	60	60000	61000	65	65000	66000
1	66	66000	67000	78	78000	79000
1	66	66000	67000	70	70000	71000
1	71	71000	72000	78	78000	79000
1	79	79000	80000	99	99000	100000
1	79	79000	80000	85	85000	86000
1	86	86000	87000	99	99000	100000
1	86	86000	87000	94	94000	95000
1	95	95000	96000	99	99000	100000
```
Each line represents one TAD of which boundary starts with the bin of the first index and ends with the bin of the last index in that line.

#
You may find more detailed introduction to SuperTAD [here](https://supertad.deepomics.org)