#!/bin/bash

set -e

INVERSE_PATH=''
gmpfrxx_dir_PATH=''
include_pkgs_PATH=''
lib_pkgs_PATH=''

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@             SMEARING EXTRACTION            @@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


# This is a bash script that allows to solve inverse problems changing
#different parameters, such as the basis, the method and the simulation
#parameters. In the following everything is commented so you can
#easily re-use the code. This code takes inputs by inputfile.txt.
#So please change it.



optimization="3"

#Choice of the method you want to use.
#If 0 -> Regular
#If 1 -> Modified
#method
boolM=1

#Choice of the basis function
#If 0  -> exp
#If 1  -> cosh
#If 2  -> omega*(cosh/sinh), used for the computation of the sphal rate
#If 3  -> ((1-exp(-beta*omega))/2) * (cosh/sinh), for the computation of the sphal rate k not zero
#basis
boolB=2

#Choice of the traget function
#If 0 Gaussian
#If 1 Pseudo-Gaussian
#target
boolTF=1

#Full covariance in the functional
#If 0 Cov(i,j)=0  if i!=j
#If 1 Cov(i,j)!=0 if i!=j
#target
boolCov=0

 
#Define Method 
if [ $boolM -eq 0 ]
then
    Method="-D BG"
fi

if [ $boolM -eq 1 ]
then
    Method="-D HLN"
fi


#Define Basis
if [ $boolB -eq 0 ]
then
    Basis="-D EXP"
fi

if [ $boolB -eq 1 ]
then
    Basis="-D COS"
fi

if [ $boolB -eq 2 ]
then
    Basis="-D COS_SPHAL"
fi

if [ $boolB -eq 3 ]
then
    Basis="-D BASIS_TEST"
fi

#Define Target function

if [ $boolTF -eq 0 ]
then
    Target="-D GAUSS"
fi

if [ $boolTF -eq 1 ]
then
    Target="-D PSEUDO_GAUSS"
fi

#Define Target function

if [ $boolCov -eq 1 ]
then
    Cov="-D COV"
fi


#Compile with all the libraries
g++ -O$optimization -std=c++14 -o main ../src/main.cpp ../lib/params.c ../lib/smear.c ../lib/statistical.c -I${gmpfrxx_dir_PATH} -L${gmpfrxx_dir_PATH} -I${include_pkgs_PATH} -L${lib_pkgs_PATH} -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -lgsl -lgslcblas $Method $Basis $Target $Cov

echo 'Ciao!'
