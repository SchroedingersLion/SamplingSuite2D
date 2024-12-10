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



class Simulation{

    public:


        Simulation()
            :  
        {

        }


        void run(input_args);



    private:

        const double stepsize;
        const double friction;
        const double temperature;
        const double N_iteration;

};





// sampler functions, implemented in .cpp file
void A_step(const double& step_size);
void B_step(const double& step_size);
void O_step(const double& a_const1, const double& a_const2, double& random_draw, std:: mt19937& RNG, std:: normal_distribution<>& normal);
void Z_step(params& parameters, const double& alpha_inv, const double& exptau);

double zeta_r;  // help var used in Sundman transform 
void Sundman_transform(const double& dtau, const double& M, const double& m, const double& r);


void ZBAOABZ(measurements& results, const double T, const double gamma, const double dtau, const double r, const double alpha, const double m, const double M, 
             const int max_iter, const int burnin, const int t_meas, const int n_dist, const int seed, const params params_init);

void BAOAB(measurements& results, const double T, const double gamma, const double dtau, const int max_iter, const int burnin, const int t_meas, const int n_dist, 
           const int seed, const params params_init);



void A_step(params& parameters, const double& step_size){
    parameters.pos_x += step_size * parameters.vel_x;
    parameters.pos_y += step_size * parameters.vel_y;
}

void B_step(params& parameters, const double& step_size){
    parameters.vel_x += step_size * parameters.force_x;
    parameters.vel_y += step_size * parameters.force_y;
}

void O_step(params& parameters, const double& a_const1, const double& a_const2, double& random_draw, std:: mt19937& RNG, std:: normal_distribution<>& normal){
    random_draw = normal(RNG);
    parameters.vel_x  =  a_const1 * parameters.vel_x  +  a_const2 * random_draw;
    random_draw = normal(RNG);
    parameters.vel_y  =  a_const1 * parameters.vel_y  +  a_const2 * random_draw;
}

void Z_step(params& parameters, const double& alpha_inv, const double& exptau){
    parameters.zeta = exptau * parameters.zeta  +  alpha_inv * (1-exptau) * (parameters.force_x * parameters.force_x  +  parameters.force_y * parameters.force_y);
}

void Sundman_transform(params& parameters, const double& dtau, const double& M, const double& m, const double& r){
    zeta_r = pow(parameters.zeta, r); 
    parameters.dt = dtau  *  m * (zeta_r + M) / (zeta_r + m);    
}




void BAOAB(measurements& results, const double T, const double gamma, const double dtau, const int max_iter, const int burnin, const int t_meas, const int n_dist, 
           const int seed, const params params_init){

    std:: cout << "Running BAOAB" << std:: endl;

    // INTEGRATOR CONSTANTS
    const double a = exp(-1*gamma*dtau);    
    const double sqrt_Ta_sq = sqrt((1-a*a)*T);
    const double dtau_half = 0.5*dtau;   

    params parameters {params_init};

    // PREPARE RNG
    std:: mt19937 twister;
    std:: seed_seq seq{1,20,3200,403,5*seed+1,12000,73667,9474+seed,19151-seed};
    std:: vector < std::uint32_t > seeds(1);
    seq.generate(seeds.begin(), seeds.end());
    twister.seed(seeds.at(0)); 

    std:: normal_distribution<> normal{0,1};
    double Rn;

    // COMPUTE INITIAL FORCES
    compute_force(parameters);

    auto t1 = std:: chrono::high_resolution_clock::now();
    std::cout<<"starting main loop"<<std::endl;

    // MAIN LOOP.
    for ( size_t i = 0;  i < max_iter;  ++i ) {

        // TAKE MEASUREMENT
        if( i % t_meas == 0  &&  i >= burnin) results.take_measurement(parameters);

        B_step(parameters, dtau_half);         

        A_step(parameters, dtau_half);        

        O_step(parameters, a, sqrt_Ta_sq, Rn, twister, normal);

        A_step(parameters, dtau_half);

        compute_force(parameters);

        B_step(parameters, dtau_half);
  
        if( i % int(1e6) == 0 ) std:: cout << "Iteration " << i << " done!" << "\n";
	
	}  // END MAIN LOOP


    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";

    return;

}