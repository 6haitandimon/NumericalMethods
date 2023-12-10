#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "gauss/gauss.h"

std::vector<std::vector<double>> get_J(double, double, double);
void Neuton(double, double, double, double, int, double);