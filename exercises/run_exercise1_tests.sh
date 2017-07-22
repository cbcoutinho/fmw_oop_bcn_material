#!/bin/bash
PATH_TO_EXERCISE1=build_exercise1/bin/exercise1

#Full matrix storage tests
echo "*************************"
echo "Full matrix storage tests"
echo "*************************"
$PATH_TO_EXERCISE1 10 10 full_matrix dir_solve
$PATH_TO_EXERCISE1 10 10 full_matrix cg_solve

#Band matrix storage tests
echo "*************************"
echo "Band matrix storage tests"
echo "*************************"
$PATH_TO_EXERCISE1 10 10 band_matrix dir_solve
$PATH_TO_EXERCISE1 10 10 band_matrix cg_solve

#Symmetric band matrix storage tests
#TO BE ADDED AS PART OF FMW OOP EXERCICE

#Sparse matrix storage tests
echo "***************************"
echo "Sparse matrix storage tests"
echo "***************************"
$PATH_TO_EXERCISE1 10 10 sparse_matrix cg_solve
