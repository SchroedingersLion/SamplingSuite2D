/*
This header file specifies the 2D problem the samplers run on as well as what observables are collected.
If you want to create a new 2D problem, you need to overwrite the force function below.
If you want to change the observables that are collected, you need to modify the measurement class below.
The code sections that require modifications for these aims are marked by the # symbol like
########
CODE
########.

Do not change anything else unless you know what you're doing.
*/

#define _USE_MATH_DEFINES
#include <stdexcept> 
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <numeric>
#include <iterator>
#include <fstream>
#include <mpi.h>



/* parameter object the samplers operate on via the "compute_force" routine below.
   it is also used by the measurement class below to obtain observable samples. */
struct params{
    double pos_x, pos_y, vel_x, vel_y, force_x{0}, force_y{0}, dt{0}, zeta{0};  // dt is the current adaptive step size for the adaptive schemes.
};



// ####### DEFINE FORCE FUNCTION ######################
// constants used in the force routine
const std:: vector <double> MB_A {-200, -100, -170, 15};  
const std:: vector <double> MB_a {-1, -1, -6.5, 0.7};
const std:: vector <double> MB_b {0, 0, 11, 0.6};
const std:: vector <double> MB_c {-10, -10, -6.5, 0.7};
const std:: vector <double> MB_x {1, 0, -0.5, -1};
const std:: vector <double> MB_y {0, 0.5, 1.5, 1};
double xdiff, ydiff, exponent;


void compute_force(params& parameters){

    parameters.force_x = parameters.force_y = 0;

    for (int i=0; i<4; ++i){
        xdiff = parameters.pos_x - MB_x[i];
        ydiff = parameters.pos_y - MB_y[i];
        exponent = MB_a[i]*xdiff*xdiff + MB_b[i]*xdiff*ydiff + MB_c[i]*ydiff*ydiff;

        parameters.force_x -= MB_A[i] * exp( exponent ) * (2*MB_a[i]*xdiff + MB_b[i]*ydiff);
        parameters.force_y -= MB_A[i] * exp( exponent ) * (2*MB_c[i]*ydiff + MB_b[i]*xdiff);
    }

}

// ####################################################




/*  This is the measurement class that is used by the samplers to obtain observables. 
    In order to modify what observables are collected, the user has to do 2 things:
    a) in the constructor, they have to specify the numbers of observables to be taken by
       adjusting the value of the variable "no_observables". 
    b) in the function "take_measurement" the user has to adjust the formulas used to compute an observable from 
       the parameters. */

class measurements {

    public:

        // constructor
        measurements(const int method_type, const int burnin, const int t_meas, const int n_dist, const int max_iter):
            method_type {method_type}, burnin {burnin}, t_meas {t_meas}, n_dist {n_dist}, max_iter {max_iter}
            {
                
                /*######## ENTER THE NUMBER OF OBSERVABLES TO COLLECT ############*/
                no_observables = 4; 
                /*################################################################*/
                

                observables.resize(no_observables);
                observable_sums.resize(no_observables);
                observable_tavgs = std:: vector <std:: vector <float>> (no_observables, std::vector <float> ((max_iter-burnin) / (n_dist*t_meas)));

                // the adaptive schemes also collect the adaptive step sizes.
                if (method_type==1){ 
                    dt_vals_raw.resize( (max_iter-burnin) / (n_dist*t_meas));
                    zeta_raw.resize( (max_iter-burnin) / (n_dist*t_meas));
                } 
            
            };



        void take_measurement(const params& parameters){

            /* ########### COMPUTE CURRENT OBSERVABLE VALUES FROM PARAMETERS ########
               The number of entries in vector "observables" must correspond to member variable "no_observables" set by the user
               in the constructor above. The reweighting for the adaptive schemes is done automatically by the "process_sample" routine. 
               The adaptive schemes will also automatically store the adaptive step size. */            
            observables[0] = parameters.pos_x;	// x-coordinate
            observables[1] = parameters.pos_y;  // y-coordinate
            observables[2] = 0.5*(parameters.vel_x*parameters.vel_x + parameters.vel_y*parameters.vel_y);           // Tkin
            observables[3] = -0.5*(parameters.pos_x*parameters.force_x + parameters.pos_y*parameters.force_y);      // Tconf
            /*########################################################################*/

            process_sample(parameters);   // computes (reweighted) time average and stores results for later print-out.

            return;

        };


        void mpi_reduction(MPI_Comm& comm, const int& rank, const int& nr_proc);  // for mpi average, implemented in .cpp file.
        void print_to_csv(const std:: string outputname);                         // for print out, implemented in .cpp file.


    private:
        std:: vector <double> observable_sums;          // sums of observable samples for time average
        std:: vector <double> observables;              // vector storing the new sample (one entry per observable)
        std:: vector <std:: vector <float>> observable_tavgs;       // stores the evolving time average for each observable
        std:: vector <std:: vector <float>> observable_printout;    // stores the mpi process average of vector "observable_tavgs" on rank 0 
        std:: vector <double> dt_vals_raw;              // stores samples of adaptive step size (only for adaptive schemes)
        std:: vector <double> zeta_raw;                 // stores samples of zeta (only for adaptive schemes)                 
        int no_observables;                             // number of observables to be taken
        const int method_type;                          // 0=constant step size scheme, 1=adaptive scheme
        int k{0}, ctr{0}, t_avg_normalizer{0};          // help variables
        const int burnin, t_meas, n_dist, max_iter;     
        double dt_sum{0};                               
        void process_sample(const params& parameters);   // add new sample to time average and store result, implemented in .cpp file.

};




// sampler functions, implemented in .cpp file
void A_step(params& parameters, const double& step_size);
void B_step(params& parameters, const double& step_size);
void O_step(params& parameters, const double& a_const1, const double& a_const2, double& random_draw, std:: mt19937& RNG, std:: normal_distribution<>& normal);
void Z_step(params& parameters, const double& alpha_inv, const double& exptau);

double zeta_r;  // help var used in Sundman transform 
void Sundman_transform(params& parameters, const double& dtau, const double& M, const double& m, const double& r);


void ZBAOABZ(measurements& results, const double T, const double gamma, const double dtau, const double r, const double alpha, const double m, const double M, 
             const int max_iter, const int burnin, const int t_meas, const int n_dist, const int seed, const params params_init);

void BAOAB(measurements& results, const double T, const double gamma, const double dtau, const int max_iter, const int burnin, const int t_meas, const int n_dist, 
           const int seed, const params params_init);
