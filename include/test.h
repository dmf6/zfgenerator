#ifndef __TEST_H
#define __TEST_H

#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "math.h"
#include <cmath>
#include <fftw3.h>
#include <string>
#include <string.h>
#include <sstream>
#include <gsl/gsl_errno.h>
/* smoothing basis spline (B-spline) */
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <boost/foreach.hpp>
#include <omp.h>
#include <vector>

#define foreach         BOOST_FOREACH
#define FMIN 0.1
#define FMAX 4.0
#define PI 3.1415926f

/* B-spline constants */
/* number of fit coefficients */
#define NCOEFFS  30
#define SPLINE_ORDER 3
    /* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
#define NBREAK   (NCOEFFS + 2 - SPLINE_ORDER)
#define NUM_POINTS 100

using namespace std;

inline int nextpow2(int x);
inline double average(double *array, int N);
inline void subtract(double *array, int N, float value);
void fit_data(ostream &os, vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d, int n);
inline double average(vector<double> *vec);
inline  void subtract(vector<double> *vec, double value);

                                  
#endif
