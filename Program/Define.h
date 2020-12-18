//#pragma once
#ifndef _DEFINE_H
#define _DEFINE_H

//using namespace std;
#include <vector>


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF PROBLEM SPECIFIC --------

//Problem specific data
std::vector <std::vector <double> > dist;	// matrix with Euclidean distance

std::vector <TNode> node;					// vector of TSP nodes




//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF BRKGA-QL --------


std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

// Input File
int debug = 1;                              // 0 - run mode      		    1 - debug mode
int ls = 1;  				                // 0 - without local search     1 - with local search
int MAXTIME = 1;                            // maximum runtime
int MAXRUNS =  1;                           // maximum number of runs of the method
unsigned MAX_THREADS = 1;            		// number of threads

// Run
char instance[50];                          // name of instance
int numLP = 0;                              // number of LP calls

// computational time (windows system)
clock_t CPUbegin,							// initial run time
        CPUbest,							// time to find the best solution in a run
        CPUend;								// end run time

// computational time (unix systems)
struct timeval Tstart, Tend, Tbest;   

//A-BRKGA
int n;                                      // size of cromossoms
int p;          	                        // size of population
double pe;              	                // fraction of population to be the elite-set
double pm;          	                    // fraction of population to be replaced by mutants
double rhoe;             	                // probability that offspring inherit an allele from elite parent

int stagnation;                             // number of generations without improvement solution
double beta;                                // perturbation intensity
double sigma;                               // pearson correlation factor

std::vector <TSol> Pop;                      	// current population
std::vector <TSol> PopInter;               		// intermediary population

TSol bestSolution;                          // best solution found in the A-BRKGA


// Sort TSol by objective function
bool sortByFitness(const TSol &lhs, const TSol &rhs) { return lhs.fo < rhs.fo; }


// Reinforcement Learning
double epsilon;                             // greed choice possibility
double lf;                                  // learning factor
double df;                                  // discount factor
double R;                                   // reward
double qTotal;                              // q*

// list of actions RL
int sizeP [] = {233, 377, 610, 987, 1597, 2584};
double Pe[] = {0.30, 0.25, 0.20, 0.15, 0.10}; 
double Pm[] = {0.25, 0.20, 0.15, 0.10, 0.05}; 
double Rhoe[] = {0.80, 0.75, 0.70, 0.65, 0.60, 0.55}; 
double Epsilon[] = {0.10, 0.15, 0.20, 0.25, 0.30};
double LF[] = {0.2, 0.4, 0.6, 0.8, 1};
double DF[] = {0.2, 0.4, 0.6, 0.8, 1};

// number of parameters in Q-table
const int par = 7;

// actions
int a0 = 0,                                 // p
    a1 = 0,                                 // pe
    a2 = 0,                                 // pm
    a3 = 0,                                 // rhoe
    a4 = 0,                                 // epsilon
    a5 = 0,                                 // lf
    a6 = 0,                                 // df
    a0Max = 0,
    a1Max = 0,
    a2Max = 0,
    a3Max = 0,
    a4Max = 0,
    a5Max = 0,
    a6Max = 0;

std::vector <std::vector <TQ> > Q;                    // Q-Table

#endif