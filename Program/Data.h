//#pragma once
#ifndef _DATA_H
#define _DATA_H

//using namespace std;
#include <vector>
#include <algorithm>    


//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------

struct TNode								// struct with node informations
{
	int id;
	double x;
	double y;
};



//------ DEFINITION OF TYPES OF BRKGA-QL --------

/***********************************************************************************
 Struct: TVecSol
 Description: struct to represent the two vector solutions (rk and problem)
************************************************************************************/
struct TVecSol
{
	int sol;                                // position of chromosome
	double rk;                              // random-key of chromosome
};


/***********************************************************************************
 Struct: TSol
 Description: struct to represent a solution problem
************************************************************************************/
struct TSol
{
    std::vector <TVecSol> vec;              // solution of the problem and random key
    double fo;                              // objetive function value
    int label;                              // defines a community solution with a number
    int similar;                            // indicates if a solution is similar to other (0 no, 1 yes)
    int flag;                               // indicates if a local search has already been performed on this solution (0 no, 1 yes)
    int promising;                          // indicates if a solution is promising to apply local search
};


/***********************************************************************************
 Struct: TQ
 Description: struct to represent a quality matrix
************************************************************************************/
struct TQ
{
    int S;
    double pVar;
    double q;
    int k;
};


#endif