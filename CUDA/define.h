#ifndef DEFINE_H
#define DEFINE_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>


const double DT = 0.025;
const int TEND = 3000;

const int NUM = 10000;

const double INP = 15.0;
const double V0 = -70;// (mV)
const double TH = -55.0;// (mV)
const double STD_N = 2.0;//standard deviation of threshold (mV)
const double TAU = 1.0;// time constant of neuron (ms)

#endif
