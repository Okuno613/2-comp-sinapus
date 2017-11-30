#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>

const double DT = 0.025;
const int TEND = 3000;

const int NUM = pow(2,20);
const int Threads = 32;

const double TAU = 1.0;// time constant of neuron (ms)

const double TH0 = 0.0;
const double a = 5.0;
const double b = 1.0;

class Simulation{
public:
    void sim();
};
#endif
