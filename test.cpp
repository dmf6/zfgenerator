#include "test.h"


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
        /* check that we are reading in current and voltage correctly */
    ofstream mysignal ("mysignal.dat");
    for (int j = 0; j < tin.size(); j++) {
        mysignal << tin[j] << " " << iin[j] << " " << vin[j] << "\n";
    }
    mysignal << endl;
    mysignal.close();
    
        /* number of elements is i...we counted as we insert into array */
    N=tin.size();
    myfile.close();

    delete [] buff;
    
    double vv = average(&vin);
    double ii = average(&iin);
   
    subtract(&vin, vv);
    subtract(&iin, ii);
    
    double dt = tin[2] - tin[1]; /* sampling frequency */
    double T = N*dt;
    double fs = N/(T/1000);
    int Nfft =  pow( 2, ceil( log( N) / log( 2 ) ) );
    
        /* USE MAX NUMBER OF THREADS IN THE SYSTEM TO PRODUCE DFT */
     fftw_plan_with_nthreads(omp_get_max_threads());
     
           /* http://www.fftw.org/doc/The-1d-Real_002ddata-DFT.html#The-1d-Real_002ddata-DFT */
        //take FFT of voltage
    fftw_complex *vout;
    vout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nfft/2+1));
    fftw_plan p1;
    p1= fftw_plan_dft_r2c_1d(N,&vin[0], vout, FFTW_ESTIMATE);
    fftw_execute(p1);

/* take FFT of current */
    fftw_complex *iout;
    iout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((Nfft/2)+1));
    fftw_plan p2;
    p2= fftw_plan_dft_r2c_1d(N,&iin[0], iout,FFTW_ESTIMATE);
    fftw_execute(p2);
    
    /* cool STL resizing of an array using vectors */
    vector<double> Is(Nfft/2+1);
    vector<double> Vs(Nfft/2+1);
   
    double *zpower, *phase, vpower, ipower, vphase, iphase;
    zpower = new double[Nfft/2+1];
    phase = new double[Nfft/2+1];
    double vre, vim, ire, iim;
    
    for (int i = 0; i < Nfft/2+1; i++) {
        vre = (vout[i][0]);
        vim = (vout[i][1]);
        vpower = (sqrt((pow(vre, 2)) + (pow(vim, 2))))/N;
        
        ire = (iout[i][0]);
        iim = (iout[i][1]);
        ipower = (sqrt((pow(ire, 2)) + (pow(iim, 2))))/N;
        Vs[i] = vpower;
        Is[i] = ipower;
        
        if (Nfft % 2 == 0) {
            if (i > 0&& i < Nfft/2) {
                Vs[i] = vpower*2;
                Is[i] = ipower*2;
            }
        }
        else {
            if (i > 0) {
                Vs[i] = vpower*2;
                Is[i] = ipower*2;
            }
        }
        
        zpower[i]  = Vs[i]/Is[i];
        vphase = atan(vim/vre);
        iphase = atan(iim/ire);
        phase[i] = vphase - iphase;
    }
    
    int fmax_idx, fmin_idx;
    fmax_idx = 0; fmin_idx=0;
    
    double *freq, *linspace;
    freq = new double[(int)((Nfft/2+1))];
    linspace = new double[(int)((Nfft/2+1))];

        /* create NFFT/2+1 numbers from 0 to 1 */
    for (int i = 0; i < (int) (Nfft/2+1); i++){
            // cout << "linspace element " << i << " is " << ((double) i)/(N/2+1) << endl;
        linspace[i] =  ((double) i)/(N/2+1);
    }
   
    for (int j = 0; j < (int)((Nfft/2+1)); j++){
            /* normalize frequency values by duration */
        freq[j] = fs/2*linspace[j];
            //cout << j << "\t" << freq[j] << endl;
        if (freq[j] <= FMIN) {
            fmin_idx = j;
        } 
        else if (freq[j] >=FMAX) {
            fmax_idx = j;
            break;
        }
    }
     
  // /* file name manipulation */
    char str1[20], str2[20], str3[20] ,str4[20];
    
    strcpy (str1, "_raw.dat");
    strcpy (str2, "_zfit.dat");
    strcpy (str3, "_zfstats.dat");
    strcpy (str4, "_phasefit.dat");
    
    char a[40], a1[40], a2[40], a3[40];
    strncpy(a, datafile, strlen(datafile)-3);
    a[strlen(datafile)-4] = '\0';    
    memcpy (a1, a, strlen(a)+1);
    memcpy (a2, a, strlen(a)+1);
    memcpy (a3, a, strlen(a)+1);
    strcat(a, str1);
    strcat(a1, str2);
    strcat(a2, str3);
    strcat(a3, str4);
    
    ofstream myfile1;
    myfile1.open (a);
        /* Only print up to 4Hz defined by FMAX */
    int n = fmax_idx;
    for(int i=0; i<n; i++) {
        myfile1 << freq[i] << "\t" << Vs[i] << "\t" << Is[i] << "\t" << zpower[i] << "\t" << phase[i] << "\n";
    }
    myfile1 << endl;
    myfile1.close();
    
     vector<double> newF (freq, freq + fmax_idx);
     vector<double> newZ (zpower, zpower + fmax_idx);
     vector<double> newP (phase, phase + fmax_idx);
     newF.erase (newF.begin(), newF.begin()+fmin_idx);
     newZ.erase (newZ.begin(), newZ.begin()+fmin_idx);
     newP.erase (newP.begin(), newP.begin()+fmin_idx);
            
     fftw_free(vout);
     fftw_free(iout);
     fftw_destroy_plan(p1);
     fftw_destroy_plan(p2);
     delete [] freq;
     delete [] linspace;
     delete [] zpower;
     delete [] phase;
     fftw_cleanup();
     fftw_cleanup_threads();


     //     compute a smoothing B-spline to get a nice impedance profile
     //    from this data set I can read all z-f stats that includes zmax, fmax, q, fwidth, z10
     ofstream myfile2, myfile3;
     myfile2.open (a1);
     myfile3.open(a3);
     double *fd = new double[NUM_POINTS];
     double *zd = new double[NUM_POINTS];
     double *pd = new double[NUM_POINTS];
     
     // vector<double> pdata (NUM_POINTS, 0);
     fit_data(myfile2, newF, newZ, fd, zd, n);
     fit_data(myfile3, newF, newP, fd, pd, n);
         
     myfile2.close();
     myfile3.close();
     vector<double> fdata (fd, fd + NUM_POINTS);
     vector<double> zdata (zd, zd + NUM_POINTS);
        /* now measure the z-f profile and write the stats to a file */
     double zmax, z0, q, z10;
     zmax = *max_element(zdata.begin(), zdata.end());
     z0 = zdata[0];
     z10 = zdata[NUM_POINTS-1];
     q = zmax - z0;
     
     int zmax_idx = distance(zdata.begin(), max_element(zdata.begin(), zdata.end()));
     double fmax = fdata[zmax_idx];
     
         // cout << "max value at " << distance(zdata.begin(), max_element(zdata.begin(), zdata.end()));
         // cout << "Resonant frequency at " << newF[zmax_idx] << endl;
     
     double  zhalf = zmax - (zmax - z0)/2;
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
     
     ofstream myfile4;
     myfile4.open (a2);
     myfile4 << zmax << "\t" << fmax<< "\t" << q << "\t" << z10 << "\t" << fhalfwidth << endl;
     myfile4.close();
     
         // cout << "Finished! Now use Grace to plot the results (http://plasma-gate.weizmann.ac.il/Grace/)" << endl;
         // cout << "raw zf data can be found in " << a << endl;
         // cout << "Smoothed zf data can be found in " << a1 << endl;
         // cout << "zf stats  (zmax, fmax, q, z10, fhalfwidth) for this profile can be found in " << a2 << endl;
            delete [] fd;
            delete [] zd;
            delete [] pd;
            
            
    return 0;
}

    /* Modify the underlying vector passed to the fit function. We would
     * like to access the fit vectors afterward so we can extract
     * statistics
     */
void fit_data(ostream &os, vector<double> &a, vector<double> &b, double *cvec, double *dvec, int n) {
    const size_t ncoeffs = NCOEFFS;
    const size_t nbreak = NBREAK;

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
    bw = gsl_bspline_alloc(SPLINE_ORDER, nbreak);
    B = gsl_vector_alloc(ncoeffs);
     
    x = gsl_vector_alloc(n);
    y = gsl_vector_alloc(n);
    X = gsl_matrix_alloc(n, ncoeffs);
    c = gsl_vector_alloc(ncoeffs);
    w = gsl_vector_alloc(n);
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    mw = gsl_multifit_linear_alloc(n, ncoeffs);
    
    double z, f, zi, sigma;
    // for (int j = 0; j < a.size(); j++) {
    //     cout << b[j] << "\n";
    // }
    
        /* this is the z-f data to be fitted */
    for (int i = 0; i < n; i++) {
        f = a[i];
        z = b[i];
        sigma = 0.1 * z;
        dy = gsl_ran_gaussian(r, sigma);     
        zi = z+dy;
        
        gsl_vector_set(x, i, f);
        gsl_vector_set(y, i, zi);
        gsl_vector_set(w, i, 1.0 / (sigma * sigma));
    }

    //     /* use uniform breakpoints on [0, 4] */
    gsl_bspline_knots_uniform(0, 35.0, bw);

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

   /* do the fit */
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

    double xi, yi, yerr;

    for (int i = 0; i <= NUM_POINTS; i++) {
            /* C++ version of linspace :) */
        xi = FMIN + i * ( (FMAX - FMIN) / NUM_POINTS );
        
        gsl_bspline_eval(xi, B, bw);
        gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
        
        os <<  xi << "\t" << yi  << "\n";
        cvec[i] = xi;
        dvec[i] = yi;
    }
    os << endl;
    

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
}


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
