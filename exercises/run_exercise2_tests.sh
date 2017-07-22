#!/bin/bash
PATH_TO_EXERCISE2=build_exercise2/bin/exercise2

#Full matrix storage tests
echo "*************************"
echo "Full matrix storage tests"
echo "*************************"
$PATH_TO_EXERCISE2 10 10 full_matrix dir_solve
$PATH_TO_EXERCISE2 10 10 full_matrix cg_solve

#Band matrix storage tests
echo "*************************"
echo "Band matrix storage tests"
echo "*************************"
$PATH_TO_EXERCISE2 10 10 band_matrix dir_solve
$PATH_TO_EXERCISE2 10 10 band_matrix cg_solve

#Symmetric band matrix storage tests
echo "***********************************"
echo "Symmetric band matrix storage tests"
echo "***********************************"
$PATH_TO_EXERCISE2 10 10 symmetric_band_matrix dir_solve
$PATH_TO_EXERCISE2 10 10 symmetric_band_matrix cg_solve


#Sparse matrix storage tests
echo "***************************"
echo "Sparse matrix storage tests"
echo "***************************"
$PATH_TO_EXERCISE2 10 10 sparse_matrix cg_solve
