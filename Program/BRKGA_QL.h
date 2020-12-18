//#pragma once
#ifndef _BRKGA_QL_H
#define _BRKGA_QL_H

//using namespace std;
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <cstring>
#include <string.h>
#include <omp.h>
#include <algorithm>
#include <sys/time.h>
#include <utility>  // pair
#include <numeric>  // iota
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip> //graph
#include <sstream> //graph
#include <fstream> //graph

#include "Data.h"
#include "Define.h"
#include "Read.h"
#include "Decoder.h"
#include "LocalSearch.h"
#include "Output.h"


//****************************** General Functions **********************************

/************************************************************************************
 Method: ABRKGA()
 Description: Apply the method A-BRKGA to solve the problem
*************************************************************************************/
void BRKGA();

/************************************************************************************
 Method: updateBestSolution()
 Description: Update the best solution found during the run
*************************************************************************************/
void updateBestSolution(TSol s);

/************************************************************************************
 Method: InitiateQTable()
 Description: Initiate the Q-Table with random values
*************************************************************************************/
void InitiateQTable();

/************************************************************************************
 Method: ChooseAction()
 Description: Choose actions and update the parameters
*************************************************************************************/
void ChooseAction(int numGeneration);

/************************************************************************************
 Method: UpdateQTable()
 Description: Update the values of Q-Table
*************************************************************************************/
void UpdateQTable();

/************************************************************************************
 Method: CREATE INITIAL SOLUTIONS
 Description: create a initial chromossom with random keys
*************************************************************************************/
TSol CreateInitialSolutions();

/************************************************************************************
 Method: PERTURBATION
 Description: perturbation similar chromossom
*************************************************************************************/
TSol Perturbation(TSol s, double beta);

/************************************************************************************
 Method: PARAMETRICUNIFORMCROSSOVER
 Description: create a new offspring with parametric uniform crossover
*************************************************************************************/
TSol ParametricUniformCrossover(int Tpe);

/************************************************************************************
 Method: PEARSON CORRELATION
 Description: calculate the Pearson correlation coefficient between two chromossoms
*************************************************************************************/
double PearsonCorrelation(std::vector <TVecSol> s1, std::vector <TVecSol> s2);

/************************************************************************************
 Metodo: IC(TSol Pop)
 Description: apply clustering method to find promising solutions in the population
*************************************************************************************/
//void IC(int tInitial, int tEnd, double sigma);
void IC(int Tpe);

/************************************************************************************
 Method: LP
 Description: Apply Label Propagation to find communities in the population
*************************************************************************************/
//void LP(std::vector < std::vector < std::pair <int, double> > > listaArestas);
void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas);

/************************************************************************************
 Method: PROMISINGLP
 Description: Find the promising solutions to represent the communities
*************************************************************************************/
//void PromisingLP(int tInitial, int tEnd);
void PromisingLP(int Tpe);

/************************************************************************************
Method: FREE MEMORY
Description: free memory of global vector
*************************************************************************************/
void FreeMemory();

//********************************** IO Functions ***********************************

/************************************************************************************
Method: WRITE RESULTS
Description: write the results in .XLS file
*************************************************************************************/
void WriteResults(double fo, double foMedia, std::vector <double> fos, float tempoMelhor, float tempoTotal);

/************************************************************************************
Method: WRITE SOLUTION
Description: write the solution in .TXT file
*************************************************************************************/
void WriteSolution(TSol s, float timeBest, float timeTotal);

/************************************************************************************
Method: WRITE SOLUTION SCREEN
Funcao: write the solution in the screen
*************************************************************************************/
void WriteSolutionScreen(TSol s, float timeBest, float timeTotal);

/************************************************************************************
 Method: RANDOMICO
 Description: Generate a double random number between min and max
*************************************************************************************/
double randomico(double min, double max);

/************************************************************************************
 Method: IRANDOMICO
 Description: Generate a int random number between min and max
*************************************************************************************/
int irandomico(int min, int max);

#endif