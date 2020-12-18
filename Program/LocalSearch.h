//#pragma once
#ifndef _LOCALSEARCH_H
#define _LOCALSEARCH_H

#include "Data.h"

/************************************************************************************
 Method: LocalSearch
 Description: RVND
*************************************************************************************/
TSol LocalSearch(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: LS1
 Description: 2-Opt
*************************************************************************************/
TSol LS1(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: LS2
 Description: NodeInsertion
*************************************************************************************/
TSol LS2(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: LS3
 Description: NodeExchange
*************************************************************************************/
TSol LS3(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: LS4
 Description: OrOpt2
*************************************************************************************/
TSol LS4(TSol s, int n, std::vector < std::vector <double> > dist);

double rand(double min, double max);
int irand(int min, int max);

#endif