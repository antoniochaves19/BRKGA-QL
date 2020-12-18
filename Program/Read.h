//#pragma once
#ifndef _READ_H
#define _READ_H

#include "Data.h"

#include <vector>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>

void ReadData(char nameTable[], int &n, std::vector <TNode>& node, std::vector <std::vector <double> >& dist);

void FreeMemoryProblem(std::vector <TNode> node, std::vector <std::vector <double> > dist);

#endif