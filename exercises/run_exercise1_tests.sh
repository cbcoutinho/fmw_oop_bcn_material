#!/bin/bash
PATH_TO_EXERCISE1=exercise1/build/bin/exercise1

#Full matrix storage tests
echo "*************************"
echo "Full matrix storage tests"
echo "*************************"
time $PATH_TO_EXERCISE1 30 30 full_matrix dir_solve
time $PATH_TO_EXERCISE1 30 30 full_matrix cg_solve

#Band matrix storage tests
echo "*************************"
echo "Band matrix storage tests"
echo "*************************"
time $PATH_TO_EXERCISE1 30 30 band_matrix dir_solve
time $PATH_TO_EXERCISE1 30 30 band_matrix cg_solve

#Symmetric band matrix storage tests
#TO BE ADDED AS PART OF FMW OOP EXERCICE

echo "*************************"
echo "Sym Band matrix storage tests"
echo "*************************"
time $PATH_TO_EXERCISE1 30 30 sym_band_matrix dir_solve
time $PATH_TO_EXERCISE1 30 30 sym_band_matrix cg_solve

#Sparse matrix storage tests
echo "***************************"
echo "Sparse matrix storage tests"
echo "***************************"
time $PATH_TO_EXERCISE1 30 30 sparse_matrix cg_solve
