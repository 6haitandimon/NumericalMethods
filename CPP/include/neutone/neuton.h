#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "gauss/gauss.h"

extern const double a;
extern const double k;

const int iter_max = 100;

double** createMatrix(int x);
void deleteMatrix(double **X, int x);
void copyMatrix(double **X, double **copyX, int x);

typedef double(*pf)(double*, double*, double, double);

double f1(double *uk1, double *uk, double t, double Tau);
double f2(double *uk1, double *uk, double t, double Tau);
double f3(double *uk1, double *uk, double t, double Tau);


double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n);
double* Neuton(double *yk_plus, double *yk, double tk, double Tau, int n);