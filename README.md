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
cd ./SuperTAD/cmake-build-release
cmake ..
make
```

## Dependencies 
* CMake 3.6+, g++, gcc
* Eigen (included in source) 

## Basic Usage  
```
usage: ./superTAD -f <input hic matrix> [options]
OPTIONS:
	-f <input path>: Input contact matrix file path
	-w <working directory path>: Working directory path (default current working directory)
	-b: Binary version
	-m: Multiple version
	-k <int>: Number of leaves in candidate coding tree (default NAN)
	-h <int>: Hierarchy number (default 2)
	--filter <true/True/TRUE/false/False/FALSE>: Filter TADs or not (default: true)
```

## Interpret output
sample output:
```
```
