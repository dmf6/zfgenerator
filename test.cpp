#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include "math.h"
#include <cmath>
#include <fstream>
#include <fftw3.h>
#include <string>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_errno.h>
/* smoothing basis spline (B-spline) */
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>


using namespace std;

#define FMIN 0.1
#define FMAX 4.0

/* B-spline constants */

    /* number of fit coefficients */
#define NCOEFFS  12
     
    /* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
#define NBREAK   (NCOEFFS - 2)


inline int nextpow2(int x);
inline double average(double *array, int N);
inline void subtract(double *array, int N, float value);

inline double average(vector<double> *vec) {
    vector<double>::iterator it;
    double sum = 0.0;
    
    for ( it=vec->begin() ; it < vec->end(); it++ ) {
        sum += *it;
     }
    return sum/(vec->size());
}

inline  void subtract(vector<double> *vec, double value) {
    for (int i = 0; i < (int)(vec->size()); i++) {
        vec->at(i) -= value;
    } 
}

/// Round up to next higher power of 2 (return x if it's already a power
/// of 2).
inline int nextpow2(int x) {
    if (x < 0)
        return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

#define PI 3.1415926f


/* program calaculates impedance data using fftw3 */
int main(int argc, char* argv[]) {
        /* User can pass in zap data file for which an impedance
         * profile is to be generated */
    char* datafile;
    if (argc < 2) {
        cout << "Usage is -in <-f filename>" << endl;
    }
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-f") {
            datafile = argv[i+1];
        }
    } 
    // cout << "Generating Z-F profile for file " << argv[2] << endl ;
    // cout << "............................................" << endl;

    int N = 0;
    
    vector <double> tin;
    vector<double> vin;
    vector<double> iin;
    
    ifstream myfile(datafile);
     
        /*buffer to hold t, v, i values read from the current line in myfile*/
    double *buff = new double[3];
    
    string line;
    if (myfile.is_open()) {
        while( getline( myfile, line ) ) {
            stringstream sstr( line );
            double value;
            int j = 0;
            while( sstr >> value) {
                buff[j] =  value;
                j++; 
            }
            tin.push_back(buff[0]);
            iin.push_back(buff[1]);
            vin.push_back(buff[2]);
            
        }
    } else {
        cout << "Cannot open file!" << endl;
    }

        /* number of elements is i...we counted as we insert into array */
    N=tin.size();
    myfile.close();
    
        //cout << N << endl;
    delete [] buff;
    
    double vv = average(&vin);
    double ii = average(&iin);
   
    subtract(&vin, vv);
    subtract(&iin, ii);
    
    double dt = tin[2] - tin[1]; /* sampling frequency */
    double T = N*dt;
    double fs = N/(T/1000);
    
    int Nfft = nextpow2(N);
    
    double *freq, *linspace;
    int fmax_idx, fmin_idx;
    fmax_idx = 0; fmin_idx=0;
    
    
    freq = new double[(int)(2*(Nfft/2+1))];
    linspace = new double[(int)(2*(Nfft/2+1))];

        /* implementation of MATLAB linspace function */
    for (int i = 0; i < (int)(2*(Nfft/2)); i++){
        linspace[i] =  (double) i/(Nfft/2);
            //cout << linspace[i] <<  "\t" << i << endl;
    }
    
    for (int j = 0; j < (int)((Nfft/2)); j++){
            /* normalize frequency values by duration */
        freq[j] = 2*(fs/2*linspace[j]);
            //cout << j << "\t" << freq[j] << endl;
        if (freq[j] <= FMIN) {
            fmin_idx = j;
        } 
        else if (freq[j] >=FMAX) {
            fmax_idx = j;
            break;
        }
    }
      
//         /* http://www.fftw.org/doc/The-1d-Real_002ddata-DFT.html#The-1d-Real_002ddata-DFT */
    /* take FFT of voltage */
    fftw_complex *vout;
    vout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*((Nfft/2)+1));
    fftw_plan p1;
    p1= fftw_plan_dft_r2c_1d(N,&vin[0],vout, FFTW_ESTIMATE);
    fftw_execute(p1);

// /* take FFT of current */
    fftw_complex *iout;
    iout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*((Nfft/2)+1));
    fftw_plan p2;
    p2= fftw_plan_dft_r2c_1d(N,&iin[0],iout,FFTW_ESTIMATE);
    fftw_execute(p2);
  
  /* file name manipulation */
    char str1[20], str2[20], str3[20];
    
    strcpy (str1, "_raw.dat");
    strcpy (str2, "_fit.dat");
    strcpy (str3, "_zfstats.dat");
    
    char a[40], a1[40], a2[40];
    strncpy(a, datafile, strlen(datafile)-3);
    a[strlen(datafile)-4] = '\0';    
    memcpy (a1, a, strlen(a)+1);
    memcpy (a2, a, strlen(a)+1);
  
     strcat(a, str1);
     strcat(a1, str2);
     strcat(a2, str3);
    
     ofstream myfile1;
     myfile1.open (a);
     double *zpower, *phase, vpower, ipower, vphase, iphase;
        /* we want to stop at 4Hz defined by FMAX */
    int n = fmax_idx;
    zpower = new double[n];
    phase = new double[n];
    for(int i=0; i<n; i++) {
        double vre = (vout[i][0]);
        double vim = (vout[i][1]);
        vpower = (sqrt((pow(vre, 2)) + (pow(vim, 2))))/N;
        double ire = (iout[i][0]);
        double iim = (iout[i][1]);
        vphase = atan(vim/vre);
        iphase = atan(iim/ire);
        phase[i] = vphase - iphase;
        ipower = (sqrt((pow(ire, 2)) + (pow(iim, 2))))/N;
        zpower[i]  = fabs((vpower)/(ipower));
        
            // write this data to a file
        myfile1 << freq[i] << "\t" << vpower << "\t" << ipower << "\t" << zpower[i] << "\t" << phase[i] << endl;
    }
     myfile1.close();

     vector<double> newF (freq, freq + fmax_idx);
     vector<double> newZ (zpower, zpower + fmax_idx);      
     newF.erase (newF.begin(),newF.begin()+fmin_idx);
     newZ.erase (newZ.begin(),newZ.begin()+fmin_idx);
    
     fftw_free(vout);
     fftw_free(iout);
     fftw_destroy_plan(p1);
     fftw_destroy_plan(p2);
     delete [] freq;
     delete [] linspace;
     delete [] zpower;
     delete [] phase;
     fftw_cleanup();
    
     // Note: A smoothing spline differs from an interpolating spline in
     //     * that the resulting curve is not required to pass through
     //     * each datapoint

     //     compute a smoothing B-spline to get a nice impedance profile
     //    from this data set I can read all z-f stats that includes zmax, fmax, q, fwidth, z10

    const size_t ncoeffs = NCOEFFS;
    const size_t nbreak = NBREAK;
        //size_t iii, jjj;
    
    gsl_bspline_workspace *bw;
    gsl_vector *B;
    double dy;
    gsl_rng *r;
    gsl_vector *x, *y;
    gsl_vector *c, *w;
    gsl_matrix *X, *cov;
    gsl_multifit_linear_workspace *mw;
    double chisq;
    
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);
    
        /* allocate a cubic bspline workspace (k = 4) */
    bw = gsl_bspline_alloc(4, nbreak);
    B = gsl_vector_alloc(ncoeffs);
     
    x = gsl_vector_alloc(n);
    y = gsl_vector_alloc(n);
    X = gsl_matrix_alloc(n, ncoeffs);
    c = gsl_vector_alloc(ncoeffs);
    w = gsl_vector_alloc(n);
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    mw = gsl_multifit_linear_alloc(n, ncoeffs);
    
    double z, f, zi, sigma;
    
        /* this is the z-f data to be fitted */
    for (int i = 0; i < n; i++) {
        f = newF[i];
        z = newZ[i];
        sigma = 0.1 * z;
        dy = gsl_ran_gaussian(r, sigma);     
        zi = z+dy;
        
        gsl_vector_set(x, i, f);
        gsl_vector_set(y, i, zi);
        gsl_vector_set(w, i, 1.0 / (sigma * sigma));
    }

    //     /* use uniform breakpoints on [0, 15] */
    gsl_bspline_knots_uniform(0.0, 20.0, bw);


        /* construct the fit matrix X */
    for (int i = 0; i < n; ++i) {
        double xi = gsl_vector_get(x, i);
        
            /* compute B_j(xi) for all j */
        gsl_bspline_eval(xi, B, bw);
        
            /* fill in row i of X */
        for (int j = 0; j < NCOEFFS; ++j) {
            double Bj = gsl_vector_get(B, j);
            gsl_matrix_set(X, i, j, Bj);
        }
    }

    //     /* do the fit */
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);


    ofstream myfile2;
    myfile2.open (a1);
    
    double xi, yi, yerr;
    const int NUM_POINTS = 100;
    vector<double> zdata (NUM_POINTS, 0);
    vector<double> fdata (NUM_POINTS, 0);
    
    for (int i = 0; i <= NUM_POINTS; i++) {
            /* C++ version of linspace :) */
        xi = FMIN + i * ( (FMAX - FMIN) / NUM_POINTS );
        
        gsl_bspline_eval(xi, B, bw);
        gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
        
        myfile2 <<  xi << "\t" << "\t" << yi  << endl;
        fdata[i] = xi;
        zdata[i] = yi;
    }  
    myfile2.close();

        /* now measure the z-f profile and write the stats to a file */
    double zmax, z0, q, z10;
    ofstream myfile3;
    myfile3.open (a2);
        //vector<int>::iterator it;
    zmax = *max_element(zdata.begin(), zdata.end());
    z0 = zdata[0];
    z10 = zdata[NUM_POINTS-1];
    q = zmax - z0;
    
    int zmax_idx = distance(zdata.begin(), max_element(zdata.begin(), zdata.end()));
    double fmax = fdata[zmax_idx];
    
    // cout << "max value at " << distance(zdata.begin(), max_element(zdata.begin(), zdata.end()));
    // cout << "Resonant frequency at " << fdata[zmax_idx] << endl;

    int zhalf = zmax - (zmax - z0)/2;
    int zhalf_idx1; int zhalf_idx2;
    zhalf_idx1 = 0;
    zhalf_idx2 = 0;
        
    vector<double>::iterator zit;
    for (zit = zdata.begin(); zit != zdata.end(); ++zit) {
        if (*zit > zhalf) {
                /* when we find the first element greater than zhalf then break */
                zhalf_idx1 = distance(zdata.begin(), zit);
            break;
        }
    }
        /* intialize iterator to end of vector and iterate backwards */
    for (zit = zdata.end(); zit != zdata.begin(); --zit) {
        if (*zit > zhalf) {
                /* when we find the first element greater than zhalf then break */
            zhalf_idx2 =distance(zdata.begin(), zit);
            break;
        }
    } 
    double fhalfwidth = fdata[zhalf_idx2] - fdata[zhalf_idx1];
    
    
    myfile3 << zmax << "\t" << fmax<< "\t" << q << "\t" << z10 << "\t" << fhalfwidth << endl;

    myfile3.close();

    // cout << "Finished! Now use Grace to plot the results (http://plasma-gate.weizmann.ac.il/Grace/)" << endl;
    // cout << "raw zf data can be found in " << a << endl;
    // cout << "Smoothed zf data can be found in " << a1 << endl;
    // cout << "zf stats  (zmax, fmax, q, z10, fhalfwidth) for this profile can be found in " << a2 << endl;

    gsl_rng_free(r);
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(w);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);

    return 0;
}







