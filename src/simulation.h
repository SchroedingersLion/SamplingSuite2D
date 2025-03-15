

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
#include "parameters.h"
// #include "measurement.h"


class Simulation{

    public:

        Simulation(const std:: string& sampler,
                   const double stepsize,
                   const double temperature,
                   const double friction,
                   const double alpha1,
                   const double alpha2,
                   const double Sundman_m,
                   const double Sundman_M,
                   const int N_iteration,
                   const int t_meas,
                   const int burnin,
                   const int rng_seed,
                   const std:: string& forcefield, 
                   const coordinate init_position, 
                   const coordinate init_velocity,
                   const double sigma_noise,
                   Measurement& results)
                :   sampler {sampler},
                    stepsize {stepsize},
                    temperature {temperature},
                    friction {friction},
                    alpha1 {alpha1},
                    alpha2 {alpha2},
                    Sundman_m {Sundman_m},
                    Sundman_M {Sundman_M},
                    N_iteration {N_iteration},
                    t_meas {t_meas},
                    burnin {burnin},
                    rng_seed {rng_seed},
                    forcefield {forcefield},
                    sigma_noise {sigma_noise},                    
                    results {results}
        {

            params.position = init_position;
            params.velocity = init_velocity;

            if (sampler == "baoab")        run_sampler = &Simulation:: run_BAOAB;
            else if (sampler == "zbaoabz") run_sampler = &Simulation:: run_ZBAOABZ;               
            else throw std:: invalid_argument( "\nInvalid sampler argument! See --help.\n" );


            // Specify 2D problem to sample:
            if (forcefield == "muellerbrown")     compute_force = &Simulation:: compute_force_MullerBrown;
            else if (forcefield == "ackley")     compute_force = &Simulation:: compute_force_Ackley;
            else if (forcefield == "rosenbrock") compute_force = &Simulation:: compute_force_Rosenbrock;
            else if (forcefield == "rosenbrock_noisy") compute_force = &Simulation:: compute_force_Rosenbrock_Noisy;
            else if (forcefield == "beale")      compute_force = &Simulation:: compute_force_Beale;
            else if (forcefield == "beale_noisy") compute_force = &Simulation:: compute_force_Beale_Noisy;
            else if (forcefield == "entropicbarrier") compute_force = &Simulation:: compute_force_EntropicBarrier;
            else if (forcefield == "star") compute_force = &Simulation:: compute_force_StarPotential;
            else if (forcefield == "harmonic_noisy") compute_force = &Simulation:: compute_force_Harmonic_Noisy;

            else throw std:: invalid_argument( "\nInvalid forcefield argument! See --help.\n" );
            

        }


        void run_MPI_simulation(int argc, char *argv[]);



    private:

        parameters params;
        const std:: string sampler;
        const double stepsize;
        const double friction;
        const double temperature;
        const double alpha1;
        const double alpha2;
        const std:: string forcefield;
        const int N_iteration;
        const int t_meas;
        const int burnin;
        Measurement& results;
        const int rng_seed;
        std:: mt19937 twister;
        const double sigma_noise;
        std:: normal_distribution<> normal{0,1};
        void (Simulation::* run_sampler)();
        void run_BAOAB();
        void run_ZBAOABZ();
        void A_step(const double stepsize);
        void B_step(const double stepsize);
        void O_step(const double a_const1, const double a_const2);
        void Z_step(const double alpha_inv, const double exptau);
        const double Sundman_M, Sundman_m;
        void Sundman_transform(const double stepsize);
        void (Simulation::* compute_force)();
        void compute_force_MullerBrown();
        void compute_force_Ackley();
        void compute_force_Rosenbrock();
        void compute_force_Rosenbrock_Noisy();
        void compute_force_Beale();
        void compute_force_Beale_Noisy();
        void compute_force_EntropicBarrier();
        void compute_force_StarPotential();
        void compute_force_Harmonic_Noisy();


};


inline void Simulation:: run_MPI_simulation(int argc, char *argv[]){

    // Set up MPI environment
    MPI_Init(&argc, &argv);				
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nr_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nr_proc);


    // Seed RNG with mpi rank
    const int seed = rank + rng_seed;

    std:: seed_seq seq{1,20,3200,403,5*seed+1,12000,73667,9474+seed,19151-seed};
    std:: vector < std::uint32_t > seeds(1);
    seq.generate(seeds.begin(), seeds.end());
    twister.seed(seeds.at(0)); 


    // Run sampler
    (this->*run_sampler)();

    results.mpi_reduction(comm, rank, nr_proc);

    // Print results.
    if (rank==0) results.print_to_csv();

    MPI_Finalize();

};





inline void Simulation:: A_step(const double step_size){
    params.position.x += step_size * params.velocity.x;
    params.position.y += step_size * params.velocity.y;
}

inline void Simulation:: B_step(const double step_size){
    params.velocity.x += step_size * params.force.x;
    params.velocity.y += step_size * params.force.y;
}


inline void Simulation:: O_step(const double a_const1, const double a_const2){
    params.velocity.x  =  a_const1 * params.velocity.x  +  a_const2 * normal(twister);
    params.velocity.y  =  a_const1 * params.velocity.y  +  a_const2 * normal(twister);
}


inline void Simulation:: Z_step(const double alpha_frac, const double exptau){
    double force_norm_sq {params.force.x * params.force.x  +  params.force.y * params.force.y};
    params.zeta = exptau * params.zeta  +  alpha_frac * (1-exptau) * force_norm_sq;
}

const double r{0.25}; 
inline void Simulation:: Sundman_transform(const double stepsize){
    double zeta_r = pow(params.zeta, r); 
    params.dt = stepsize  *  Sundman_m * (zeta_r + Sundman_M) / (zeta_r + Sundman_m);    
}



inline void Simulation:: run_BAOAB(){

    std:: cout << "Running BAOAB" << std:: endl;

    // INTEGRATOR CONSTANTS
    const double a = exp(-1*friction*stepsize);    
    const double sqrt_Ta_sq = sqrt((1-a*a)*temperature);
    const double stepsize_half = 0.5*stepsize;   

    // COMPUTE INITIAL FORCES
    (this->*compute_force)();

    auto t1 = std:: chrono::high_resolution_clock::now();
    std::cout<<"starting main loop"<<std::endl;

    // MAIN LOOP.
    for ( size_t i = 0;  i <= N_iteration;  ++i ) {

        // TAKE MEASUREMENT
        if( i % t_meas == 0  &&  i >= burnin) results.take_measurement(params);

        B_step(stepsize_half);         

        A_step(stepsize_half);        

        O_step(a, sqrt_Ta_sq);

        A_step(stepsize_half);

        (this->*compute_force)();

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
    (this->*compute_force)();

    // INIT ZETA
    force_norm_sq = params.force.x*params.force.x + params.force.y*params.force.y;
    params.zeta = force_norm_sq; 

    // OBTAIN FIRST ADAPTIVE STEP SIZE
    Sundman_transform(stepsize);
    dt_half = 0.5*params.dt;

    auto t1 = std:: chrono::high_resolution_clock::now();
    std::cout<<"starting main loop"<<std::endl;

    // MAIN LOOP
    for ( size_t i = 0;  i <= N_iteration;  ++i ) {
        
        // TAKE MEASUREMENT
        if( i % t_meas == 0  &&  i >= burnin) results.take_measurement(params);

        Z_step(alpha_frac, exptau_half);

        Sundman_transform(stepsize);    // update step size
        dt_half = 0.5*params.dt;

        B_step(dt_half);

        A_step(dt_half);

        a = pow(e_min_gamma, params.dt);       // update constants for O-step with new dt
        sqrt_Ta_sq = sqrt((1-a*a)*temperature);
        O_step(a, sqrt_Ta_sq);

        A_step(dt_half);

        (this->*compute_force)();

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

// MullerBrown
const std:: vector <double> MullerBrown_A {-200, -100, -170, 15};  
const std:: vector <double> MullerBrown_a {-1, -1, -6.5, 0.7};
const std:: vector <double> MullerBrown_b {0, 0, 11, 0.6};
const std:: vector <double> MullerBrown_c {-10, -10, -6.5, 0.7};
const std:: vector <double> MullerBrown_x {1, 0, -0.5, -1};
const std:: vector <double> MullerBrown_y {0, 0.5, 1.5, 1};
double MullerBrown_xdiff, MullerBrown_ydiff, MullerBrown_exponent;

inline void Simulation:: compute_force_MullerBrown(){

    params.force.x = params.force.y = 0;  // Reset forces.

    for (int i=0; i<4; ++i){
        MullerBrown_xdiff = params.position.x - MullerBrown_x[i];
        MullerBrown_ydiff = params.position.y - MullerBrown_y[i];
        MullerBrown_exponent =   MullerBrown_a[i]*MullerBrown_xdiff*MullerBrown_xdiff 
                               + MullerBrown_b[i]*MullerBrown_xdiff*MullerBrown_ydiff 
                               + MullerBrown_c[i]*MullerBrown_ydiff*MullerBrown_ydiff;

        params.force.x -= MullerBrown_A[i]*exp( MullerBrown_exponent ) * (2*MullerBrown_a[i]*MullerBrown_xdiff + MullerBrown_b[i]*MullerBrown_ydiff);
        params.force.y -= MullerBrown_A[i]*exp( MullerBrown_exponent ) * (2*MullerBrown_c[i]*MullerBrown_ydiff + MullerBrown_b[i]*MullerBrown_xdiff);
    }

}


// Ackley
const double Ackley_2pi {2*M_PI};
double Ackley_distconst;
double Ackley_x, Ackley_y;

inline void Simulation:: compute_force_Ackley(){

    // Fill help constants.
    Ackley_x = params.position.x;
    Ackley_y = params.position.y;
    Ackley_distconst = sqrt(0.5*(Ackley_x*Ackley_x + Ackley_y*Ackley_y));
    
    // Compute force.
    params.force.x = -2*Ackley_x*exp(-0.2*Ackley_distconst)/Ackley_distconst 
                      - M_PI*sin(Ackley_2pi*Ackley_x)*exp(0.5*(cos(Ackley_2pi*Ackley_x)+cos(Ackley_2pi*Ackley_y)));
    
    params.force.y = -2*Ackley_y*exp(-0.2*Ackley_distconst)/Ackley_distconst 
                      - M_PI*sin(Ackley_2pi*Ackley_y)*exp(0.5*(cos(Ackley_2pi*Ackley_x)+cos(Ackley_2pi*Ackley_y)));

}


// Rosenbrock
inline void Simulation:: compute_force_Rosenbrock(){

    params.force.x = 400*params.position.x*(params.position.y - params.position.x*params.position.x)
                     + 2*(1-params.position.x);
    
    params.force.y = -200*(params.position.y - params.position.x*params.position.x);

}



// Beale
// double Beale_var1, Beale_var2, Beale_var3;
// double y_sq;

// inline void Simulation:: compute_force_Beale(){

//     y_sq = params.position.y*params.position.y;
//     Beale_var1 = 2*(1.5   - params.position.x + params.position.x*params.position.y);
//     Beale_var2 = 2*(2.25  - params.position.x + params.position.x*y_sq);
//     Beale_var3 = 2*(2.625 - params.position.x + params.position.x*y_sq*params.position.y);

//     params.force.x = -Beale_var1*(-1+params.position.y) - Beale_var2*(-1+y_sq) - Beale_var3*(-1+y_sq*params.position.y);
//     params.force.y = -params.position.x * (Beale_var1 + 2*params.position.y*Beale_var2 + 3*y_sq*Beale_var3);

// }


// Beale confined
double Beale_var1, Beale_var2, Beale_var3;
double y_sq, x_to_5, y_to_5, Beale_confined_exponential;
constexpr double Beale_confined_a {0.00001}, Beale_confined_b {0.3};
inline void Simulation:: compute_force_Beale(){

    y_sq = params.position.y*params.position.y;
    
    x_to_5 = params.position.x * params.position.x * params.position.x * params.position.x * params.position.x;
    y_to_5 = y_sq * y_sq *params.position.y;
    Beale_confined_exponential = exp( Beale_confined_a * (x_to_5*params.position.x + y_to_5*params.position.y) );
    
    Beale_var1 = 2*(1.5   - params.position.x + params.position.x*params.position.y);
    Beale_var2 = 2*(2.25  - params.position.x + params.position.x*y_sq);
    Beale_var3 = 2*(2.625 - params.position.x + params.position.x*y_sq*params.position.y);

    params.force.x = -Beale_var1*(-1+params.position.y) - Beale_var2*(-1+y_sq) - Beale_var3*(-1+y_sq*params.position.y) - 6*Beale_confined_a*Beale_confined_b*x_to_5*Beale_confined_exponential;
    params.force.y = -params.position.x * (Beale_var1 + 2*params.position.y*Beale_var2 + 3*y_sq*Beale_var3) - 6*Beale_confined_a*Beale_confined_b*y_to_5*Beale_confined_exponential;

}



// Noisy Beale
const size_t noise_freq {500000};
size_t noise_ctr {1};
bool noise_switch {0};
inline void Simulation:: compute_force_Beale_Noisy(){

    y_sq = params.position.y*params.position.y;
    
    x_to_5 = params.position.x * params.position.x * params.position.x * params.position.x * params.position.x;
    y_to_5 = y_sq * y_sq *params.position.y;
    Beale_confined_exponential = exp( Beale_confined_a * (x_to_5*params.position.x + y_to_5*params.position.y) );

    Beale_var1 = 2*(1.5   - params.position.x + params.position.x*params.position.y);
    Beale_var2 = 2*(2.25  - params.position.x + params.position.x*y_sq);
    Beale_var3 = 2*(2.625 - params.position.x + params.position.x*y_sq*params.position.y);

    params.force.x = -Beale_var1*(-1+params.position.y) - Beale_var2*(-1+y_sq) - Beale_var3*(-1+y_sq*params.position.y) - 6*Beale_confined_a*Beale_confined_b*x_to_5*Beale_confined_exponential;
    params.force.y = -params.position.x * (Beale_var1 + 2*params.position.y*Beale_var2 + 3*y_sq*Beale_var3) - 6*Beale_confined_a*Beale_confined_b*y_to_5*Beale_confined_exponential;

    // if (noise_switch){
        params.force.x += sigma_noise*normal(twister);
        params.force.y += sigma_noise*normal(twister);
    // }
    ++noise_ctr;
    if (noise_ctr%noise_freq==0) noise_switch = !noise_switch;

}




// Entropic Channel
double EntropicBarrier_x2;
double EntropicBarrier_denominator;

inline void Simulation:: compute_force_EntropicBarrier(){

    EntropicBarrier_x2 = params.position.x * params.position.x;
    EntropicBarrier_denominator = 1+10*EntropicBarrier_x2*EntropicBarrier_x2;

    params.force.x = 4000 * (EntropicBarrier_x2*params.position.x*params.position.y*params.position.y) / (EntropicBarrier_denominator*EntropicBarrier_denominator) - 0.004*(EntropicBarrier_x2-9)*params.position.x;

    params.force.y = -200 * params.position.y / EntropicBarrier_denominator;

}



// Noisy Rosenbrock
constexpr double sigma_rosenbrock {2};
inline void Simulation:: compute_force_Rosenbrock_Noisy(){

    params.force.x = 400*params.position.x*(params.position.y - params.position.x*params.position.x)
                     + 2*(1-params.position.x)
                     + sigma_rosenbrock*normal(twister);
    
    params.force.y = -200*(params.position.y - params.position.x*params.position.x) + sigma_rosenbrock*normal(twister);

}



inline void Simulation:: compute_force_StarPotential(){

    params.force.x = -2*params.position.x - 2000*params.position.x*params.position.y*params.position.y;
    params.force.y = -2*(1+1000*params.position.x*params.position.x) * params.position.y;

}

constexpr double wx_harmonic = 1;
constexpr double wy_harmonic = 100;
inline void Simulation:: compute_force_Harmonic_Noisy(){

    params.force.x = -2*wx_harmonic*params.position.x;
    params.force.y = -2*wy_harmonic*params.position.y - params.position.y*sigma_noise*abs(normal(twister));

}