//#pragma once
#ifndef _DECODER_H
#define _DECODER_H

#define INFINITO 999999999

#include <math.h>
#include "Data.h"


/************************************************************************************
 Method: Decoder()
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
TSol Decoder(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: Dec1
 Description: standard decoder 
*************************************************************************************/
TSol Dec1(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: Dec2
 Description: standard decoder 
*************************************************************************************/
TSol Dec2(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: Dec3
 Description: standard decoder 
*************************************************************************************/
TSol Dec3(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: Dec4
 Description: standard decoder 
*************************************************************************************/
TSol Dec4(TSol s, int n, std::vector < std::vector <double> > dist);

/************************************************************************************
 Method: Dec5
 Description: standard decoder 
*************************************************************************************/
TSol Dec5(TSol s, int n, std::vector < std::vector <double> > dist);

#endif