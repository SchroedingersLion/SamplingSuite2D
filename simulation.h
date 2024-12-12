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

struct coordinate{
    double x{0}, y{0};
}

struct parameters{
    coordinate position, velocity, force;
    double zeta, dt;    // used only by SamAdams.
}


class Simulation{

    public:


        Simulation(const std:: string& sampler,
                   const double stepsize,
                   const double temperature,
                   const double friction,
                   const double alpha1,
                   const double alpha2,
                   const int N_iteration,
                   const int t_meas,
                   const int burnin,
                   const std:: string& forcefield, 
                   const coordinate init_position, 
                   const coordinate init_velocity,
                   const bool time_average,
                   measurement& results)
                :   sampler {sampler},
                    stepsize {stepsize},
                    temperature {temperature},
                    friction {friction},
                    alpha1 {alpha1},
                    alpha2 {alpha2},
                    N_iteration {N_iteration},
                    t_meas {t_meas},
                    burnin {burnin},
                    forcefield {forcefield},
                    time_average {time_average},
                    results {results}
        {

            parameters.position = init_position;
            parameters.velocity = init_velocity;

            if (sampler="baoab")        run_sampler = &(Simulation:: run_BAOAB);
            else if (sampler="zbaoabz") run_sampler = &(Simulation:: run_ZBAOABZ);               
            else throw std:: invalid_argument( "\nInvalid sampler argument! See --help.\n" );


            // Specify 2D problem to sample:
            switch (forcefield){
                
                case "MullerBrown":
                    compute_force = &(Simulation:: compute_force_MullerBrown);
                    break;
                
                default:
                  throw std:: invalid_argument( "\nInvalid forcefield argument! See --help.\n" );
                  break;
            } 
            

        }


        void run_MPI_simulation(int argc, char *argv[]);



    private:

        parameters parameters;
        const double stepsize;
        const double friction;
        const double temperature;
        const double N_iteration;
        const int t_meas;
        const int burnin;
        const bool time_average;
        measurement& results;
        std:: mt19937 twister;
        std:: normal_distribution<> normal{0,1};
        void (Simulation::* run_sampler)();
        void run_BAOAB();
        void run_ZBAOABZ();
        void A_step(const double stepsize);
        void B_step(const double stepsize);
        void O_step(const double a_const1, const double a_const2);
        void Z_step(const double alpha_inv, const double exptau);
        void Sundman_transform(const double stepsize);
        void (Simulation::* compute_force)();
        void compute_force_MullerBrown();

};


inline void Simulation:: run_MPI_simulation(int argc, char *argv[]){

    // Set up MPI environment
    MPI_Init(&argc, &argv);				
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nr_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nr_proc);


    // Seed RNG with mpi rank
    const int seed = rank;

    std:: seed_seq seq{1,20,3200,403,5*seed+1,12000,73667,9474+seed,19151-seed};
    std:: vector < std::uint32_t > seeds(1);
    seq.generate(seeds.begin(), seeds.end());
    twister.seed(seeds.at(0)); 


    // Run sampler
    run_sampler();

    results.mpi_reduction(comm, rank, nr_proc);

    MPI_Finalize();

};





inline void Simulation:: A_step(const double step_size){
    parameters.position.x += step_size * parameters.velocity.x;
    parameters.position.y += step_size * parameters.velocity.y;
}

inline void Simulation:: B_step(const double step_size){
    parameters.velocity.x += step_size * parameters.force.x;
    parameters.velocity.y += step_size * parameters.force.y;
}


inline void Simulation:: O_step(const double a_const1, const double a_const2){
    parameters.velocity.x  =  a_const1 * parameters.velocity.x  +  a_const2 * normal(RNG);
    parameters.velocity.y  =  a_const1 * parameters.velocity.y  +  a_const2 * normal(RNG);
}


inline void Simulation:: Z_step(const double alpha_frac, const double exptau){
    double force_norm_sq {parameters.force.x * parameters.force.x  +  parameters.force.y * parameters.force.y};
    parameters.zeta = exptau * parameters.zeta  +  alpha_frac * (1-exptau) * force_norm_sq;
}

const double M{10}, m{0.1}, r{0.25}; // leave them here for convenience.
inline void Simulation:: Sundman_transform(const double stepsize){
    double zeta_r = pow(parameters.zeta, r); 
    parameters.dt = stepsize  *  m * (zeta_r + M) / (zeta_r + m);    
}



inline void Simulation:: run_BAOAB(){

    std:: cout << "Running BAOAB" << std:: endl;

    // INTEGRATOR CONSTANTS
    const double a = exp(-1*friction*stepsize);    
    const double sqrt_Ta_sq = sqrt((1-a*a)*temperature);
    const double stepsize_half = 0.5*stepsize;   

    // COMPUTE INITIAL FORCES
    compute_force();

    auto t1 = std:: chrono::high_resolution_clock::now();
    std::cout<<"starting main loop"<<std::endl;

    // MAIN LOOP.
    for ( size_t i = 0;  i < N_iterations;  ++i ) {

        // TAKE MEASUREMENT
        if( i % t_meas == 0  &&  i >= burnin) results.take_measurement(parameters);

        B_step(stepsize_half);         

        A_step(stepsize_half);        

        O_step(a, sqrt_Ta_sq);

        A_step(stepsize_half);

        compute_force();

        B_step(stepsize_half);
  
        if( i % int(1e6) == 0 ) std:: cout << "Iteration " << i << " done!" << "\n";
	
	}  // END MAIN LOOP


    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";

    return;

}



inline void Simulation:: run_ZBAOABZ(){
    
    std:: cout << "Running ZBAOABZ" << std:: endl;

    // INTEGRATOR CONSTANTS
    const double e_min_gamma {exp(-friction)};
    const double exptau_half{exp(-alpha1*0.5*stepsize)};
    double sqrt_Ta_sq;
    double a;

    const double alpha_frac {alpha2/alpha1};
    double dt_half;
    double force_norm_sq;

    // COMPUTE INITIAL FORCES
    compute_force();

    // INIT ZETA
    force_norm_sq = parameters.force.x*parameters.force.x + parameters.force.y*parameters.force.y;
    parameters.zeta = force_norm_sq; 

    // OBTAIN FIRST ADAPTIVE STEP SIZE
    Sundman_transform(stepsize);
    dt_half = 0.5*parameters.dt;

    auto t1 = std:: chrono::high_resolution_clock::now();
    std::cout<<"starting main loop"<<std::endl;

    // MAIN LOOP
    for ( size_t i = 0;  i < N_iteration;  ++i ) {
        
        // TAKE MEASUREMENT
        if( i % t_meas == 0  &&  i >= burnin) results.take_measurement(parameters);

        Z_step(alpha_frac, exptau_half);

        Sundman_transform(stepsize);    // update step size
        dt_half = 0.5*parameters.dt;

        B_step(dt_half);

        A_step(dt_half);

        a = pow(e_min_gamma, parameters.dt);       // update constants for O-step with new dt
        sqrt_Ta_sq = sqrt((1-a*a)*T);
        O_step(a, sqrt_Ta_sq);

        A_step(dt_half);

        compute_force();

        B_step(dt_half);
  
        Z_step(alpha_frac, exptau_half);

        Sundman_transform(stepsize);     // update step size
		
        if( i % int(1e6) == 0 ) std:: cout << "Iteration " << i << " done!" << "\n";
	
    }  // END MAIN LOOP


    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";

    return;
}



// FORCES

// constants used in the force routine
const std:: vector <double> MullerBrown_A {-200, -100, -170, 15};  
const std:: vector <double> MullerBrown_a {-1, -1, -6.5, 0.7};
const std:: vector <double> MullerBrown_b {0, 0, 11, 0.6};
const std:: vector <double> MullerBrown_c {-10, -10, -6.5, 0.7};
const std:: vector <double> MullerBrown_x {1, 0, -0.5, -1};
const std:: vector <double> MullerBrown_y {0, 0.5, 1.5, 1};
double MullerBrown_xdiff, MullerBrown_ydiff, MullerBrown_exponent;

inline void Simulation:: compute_force_MullerBrown(){

    parameters.force.x = parameters.force.y = 0;  // Reset forces.

    for (int i=0; i<4; ++i){
        MullerBrown_xdiff = parameters.position.x - MullerBrown_x[i];
        MullerBrown_ydiff = parameters.position.y - MullerBrown_y[i];
        MullerBrown_exponent =   MullerBrown_a[i]*MullerBrown_xdiff*MullerBrown_xdiff 
                               + MullerBrown_b[i]*MullerBrown_xdiff*MullerBrown_ydiff 
                               + MullerBrown_c[i]*MullerBrown_ydiff*MullerBrown_ydiff;

        parameters.force.x -= MullerBrown_A[i]*exp( MullerBrown_exponent ) * (2*MullerBrown_a[i]*MullerBrown_xdiff + MullerBrown_b[i]*MullerBrown_ydiff);
        parameters.force.y -= MullerBrown_A[i]*exp( MullerBrown_exponent ) * (2*MullerBrown_c[i]*MullerBrown_ydiff + MullerBrown_b[i]*MullerBrown_xdiff);
    }

}


