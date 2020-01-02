# SuperTAD

## Install  
use git:  
```
git clone git@gitlab.deepomics.org:mengbo/SuperTAD_CPP.git SuperTAD
```
or download from source
```
wget https://supertad.deepomics.org/home/download_src SuperTAD.tar.gz
tar -xzvf SuperTAD.tar.gz
```
then
```
cd ./SuperTAD
cmake .
make
```

## Dependencies 
* CMake 3.5.1+, g++, gcc
* Eigen (included in the source)

## Basic Usage  
```
usage: ./superTAD -f <input hic matrix> [options]
OPTIONS:
	-f <input path>: Input contact matrix file path
	-b: Binary version (default)
	-m: Multiple version
	-k <int>: Number of leaves in candidate coding tree (default NAN)
	-h <int>: Hierarchy number (default 2)
	--filter <true/True/TRUE/false/False/FALSE>: Filter TADs or not (default: true)
    -v/--verbose: Print verbose\n";
```

## Interpret output
sample output:
```
135 136 137 138 139 140 141 142 
128 129 130 131 132 133 134 
112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 
90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 
80 81 82 83 84 85 86 87 88 89 
65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 
55 56 57 58 59 60 61 62 63 64 
46 47 48 49 50 51 52 53 54 
38 39 40 41 42 43 44 45 
26 27 28 29 30 31 32 33 34 35 36 37 
11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
1 2 3 4 5 6 7 8 9 10
```
each line represents one TAD of which boundary starts with the bin of the first index and ends with the bin of the last index in that line.

You may access SuperTAD at [supertad.deepomics.org](https://supertad.deepomics.org)