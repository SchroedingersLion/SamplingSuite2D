/*
C++ source file to sample 2-dimensional potential energy surfaces specified by a force function.
The available samplers are BAOAB, ZBAOABZ, and OABZBAO. 

The code needs to be compiled with mpi wrappers. Using OpenMPI and the Gnu compiler, the command is
"mpicxx -O3 -o adam_sampling_2D.exe adam_sampling_2D.cpp".

The .exe needs to be run with 4 command line arguments: The first is an integer specifying the method, the second specifies the
step size parameter dtau, the third specifies friction gamma, the fourth zeta-friction alpha (ignored in BAOAB).
Example: "mpirun -n 10 sampling_2D.exe 2 0.01 1 10" to average over 10 independent trajectories of ZBAOABZ with dtau=0.01,
gamma=1, alpha=10.

The execution generates an output file of N+1 columns (BAOAB) or N+2 columns (adaptive scheme),
where N is the number of observables to be taken.
1st column:      iteration count
2nd column:      averaged first observable
...
(N+1)-th column: averaged N-th observable 
(N+2)-th column: samples of adaptive step size of mpi rank 0 (only for the adaptive schemes)

In order to change the sampling parameters apart from what can be specified by the command line arguments,
the user needs to modify the main() function below.
In order to change the 2D problem the samplers run on or change the observables that are to be taken,
the user needs to modify the header file "adam_sampling_2D.h" (see explanation on top of that file).

All sections the user is encouraged to modify for these aims, in this file as well as the header, are marked by # like
#####
// code
##### 

*/


#include "adam_sampling_2D.h"   // specifies force function as well as which obervables to take



// MAIN FUNCTION
int main (int argc, char *argv[]){

    // input arguments
    const int method = atoi(argv[1]);       // 0=BAOAB, 1=OABZBAO, 2=ZBAOABZ,  
    const double dtau = atof(argv[2]);      // step size parameter dtau  
    const double gamma = atof(argv[3]);
    const double alpha = atof(argv[4]);


    // ############ PARAMETERS TO SET ################################

    const double T = 1;         // temperature

    // parameters for Adam samplers (if used)
    const double r = 0.25;
    const double m = 0.1;
    const double M = 10;

    int max_iter = 5e7;     // iteration number

    const int burnin = 0;    // discard first burnin samples
    const int t_meas = 2;    // take observable sample and add it to moving average any t_meas iterations
    const int n_dist = 250;  // store and print any n_dist taken sample

    // initial conditions
    const params params_init{-0.6, 1.4, 0, 0};  // {pos x, pos y, vel x, vel y}

    // #################################################################


    // measurement object the samplers will operate on
    const int method_type {method==0 ? 0 : 1};
    measurements results(method_type, burnin, t_meas, n_dist, max_iter);

    // output file name
    std:: string outputfile;
    if (method == 0)        outputfile = std::string("BAOAB_dtau")+argv[2]+ "_gam" + argv[3] + ".csv";
    else if (method == 1)   outputfile = std::string("ZBAOABZ_dtau")+argv[2]+ "_gam" + argv[3] + "_alpha" + argv[4] + ".csv";

    // INIT MPI
    MPI_Init(&argc, &argv);				
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nr_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nr_proc);

    const int seed = rank;

    // RUN SAMPLER
    if (method == 0)      BAOAB(results, T, gamma, dtau, max_iter, burnin, t_meas, n_dist, seed, params_init);
    else if (method == 1) ZBAOABZ(results, T, gamma, dtau, r, alpha, m, M, max_iter, burnin, t_meas, n_dist, seed, params_init);
    else                  throw std::invalid_argument( "received unassigned method number." );

    // MPI AVERAGE
    results.mpi_reduction(comm, rank, nr_proc);

    // PRINT TO FILE
    if( rank==0 ) results.print_to_csv(outputfile);        


    MPI_Finalize();
    return 0;

}





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



void ZBAOABZ(measurements& results, const double T, const double gamma, const double dtau, const double r, const double alpha, const double m, const double M, 
             const int max_iter, const int burnin, const int t_meas, const int n_dist, const int seed, const params params_init){

    std:: cout << "Running ZBAOABZ" << std:: endl;

    // INTEGRATOR CONSTANTS
    const double e_min_gamma {exp(-gamma)};
    const double exptau_half{exp(-alpha*0.5*dtau)};
    double sqrt_Ta_sq;
    double a;

    const double alpha_inv {1/alpha};
    double dt_half;
    double force_norm_sq;
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

    // INIT ZETA
    force_norm_sq = parameters.force_x*parameters.force_x + parameters.force_y*parameters.force_y;
    parameters.zeta = force_norm_sq; 

    // OBTAIN FIRST ADAPTIVE STEP SIZE
    Sundman_transform(parameters, dtau, M, m, r);
    dt_half = 0.5*parameters.dt;

    auto t1 = std:: chrono::high_resolution_clock::now();
    std::cout<<"starting main loop"<<std::endl;

    // MAIN LOOP
    for ( size_t i = 0;  i < max_iter;  ++i ) {
        
        // TAKE MEASUREMENT
        if( i % t_meas == 0  &&  i >= burnin) results.take_measurement(parameters);

        Z_step(parameters, alpha_inv, exptau_half);

        Sundman_transform(parameters, dtau, M, m, r);;    // update step size
        dt_half = 0.5*parameters.dt;

        B_step(parameters, dt_half);

        A_step(parameters, dt_half);

        a = pow(e_min_gamma, parameters.dt);       // update constants for O-step with new dt
        sqrt_Ta_sq = sqrt((1-a*a)*T);
        
        O_step(parameters, a, sqrt_Ta_sq, Rn, twister, normal);

        A_step(parameters, dt_half);

        compute_force(parameters);

        B_step(parameters, dt_half);
  
        Z_step(parameters, alpha_inv, exptau_half);

        Sundman_transform(parameters, dtau, M, m, r);     // update step size
		
        if( i % int(1e6) == 0 ) std:: cout << "Iteration " << i << " done!" << "\n";
	
    }  // END MAIN LOOP


    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";

    return;
}





// implementation of measurement class routines that the user does not need to modify.

void measurements:: print_to_csv(const std:: string outputname){
    
    std:: ofstream file{outputname};
    std:: cout << "Writing to file...\n";

    // annoying logic to obtain first iteration index at which observable will be printed to file 
    // (depends on variables burnin, t_meas, and n_dist).
    int first_index;  
    if (burnin >= t_meas) first_index = (burnin % t_meas) == 0  ?  burnin :  burnin - (burnin % t_meas) + t_meas;
    else if ( (0<burnin) && (burnin<t_meas) ) first_index = t_meas;
    else if ( burnin == 0 ) first_index = 0;
    else throw std::invalid_argument( "Something wrong in write function." );

    // print results to file
    for ( size_t i = 0; i<observable_printout[0].size(); ++i )
    {
        file << first_index + i*t_meas*n_dist << " ";
        for ( size_t j = 0; j<no_observables; ++j ){
            file << observable_printout[j][i] << " ";  
        }
        
        if (dt_vals_raw.size() > 0) file << dt_vals_raw[i] << " ";
        if (zeta_raw.size() > 0) file << zeta_raw[i] << " ";
        
        file << "\n";
    }

    file.close();

    return;

}


void measurements:: mpi_reduction(MPI_Comm& comm, const int& rank, const int& nr_proc){   // MPI AVERAGE ROUTINE
    
    // resize output array to store averaged results in
    observable_printout.resize(no_observables);
    const int row_size = observable_tavgs[0].size();
    if (rank==0){
        for (size_t i=0; i<no_observables; ++i){
            observable_printout[i].resize( row_size ); 
        }    
    }

    // collect results from processes
    for (size_t i=0; i<no_observables; ++i) MPI_Reduce( &observable_tavgs[i][0], &observable_printout[i][0], row_size, MPI_FLOAT, MPI_SUM, 0, comm);  

    if( rank==0 ){
        for (size_t i=0; i<no_observables; ++i){
            for (size_t j=0;  j<row_size; ++j){
                observable_printout[i][j] /= nr_proc;     // divide by no. of processes to obtain averages
            }
        }
    }

    return;      
};


void measurements:: process_sample(const params& parameters){

    if (method_type==0){   // constant step size scheme, eg. BAOAB.

        for (int i=0; i<no_observables; ++i) observable_sums[i] += observables[i];  // add to sum for time average

        ++t_avg_normalizer;

        // store current moving average value for print-out, but only every n_dist steps        
        if (ctr % n_dist == 0){                                             
            for (int i=0; i<no_observables; ++i) observable_tavgs[i][k] = observable_sums[i] / t_avg_normalizer;
            ++k;
        }          
    }

    else if (method_type==1){   // adaptive step size scheme.

        for (int i=0; i<no_observables; ++i) observable_sums[i] += observables[i] * parameters.dt;  // reweighting
        dt_sum += parameters.dt;

        if (ctr % n_dist == 0){                                                 
            for (int i=0; i<no_observables; ++i) observable_tavgs[i][k] = observable_sums[i] / dt_sum;
            dt_vals_raw[k] = parameters.dt;
            zeta_raw[k] = parameters.zeta;
            ++k;
        }         
    }

    ++ctr;

}
