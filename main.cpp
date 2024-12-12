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


#include "simulation.h"   // specifies force function as well as which obervables to take
#include "measurement.h"




int main(int argc, char *argv[]){

    // Process command line arguments.
    auto [parse_result, options] = parseCommandLine(argc, argv); // Parse command line.
    if (parse_result.count("help")) {std:: cout << options.help() << std:: endl; return 0;}
    ParsedValues vals = processParsedValues(parse_result);

    // Set up measurement object.
    bool method_type = vals.sampler=="ZBAOABZ" ? 1 : 0;
    results = measurements(method_type, vals.burnin, vals.t_meas, vals.n_dist, vals.max_iter);
    
    // Set up simulation object.
    simulation simu(vals.sampler, 
                    vals.stepsize, 
                    vals.temperature, 
                    vals.friction, 
                    vals.alpha1, 
                    vals.alpha2, 
                    vals.N_iteration, 
                    vals.t_meas, 
                    vals.burnin, 
                    vals.forcefield, 
                    vals.init_position, 
                    vals.init_velocity,
                    vals.time_average,
                    results);
    
    // Run simulation.
    simu.run_MPI_simulation();

    // Print results.
    results.print_results(vals.output_name);


    return 0;

}


// MAIN FUNCTION
int main (int argc, char *argv[]){

    // input arguments



    // ############ PARAMETERS TO SET ################################

    const double T = 1;         // temperature

    // parameters for Adam samplers (if used)
    const double r = 0.25;
    const double m = 0.1;
    const double M = 10;

    int max_iter = 5e7;     // iteration number

    const int burnin = 0;   // discard first burnin samples
    const int t_meas = 2;   // take observable sample and add it to moving average any t_meas iterations
    const int n_dist = 250; // store and print any n_dist taken sample

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




    





