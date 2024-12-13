/* THIS DESCRIPTION NEEDS TO BE UPDATED
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

#include "parameters.h"
#include "measurement.h"
#include "simulation.h"   // specifies force function as well as which obervables to take
#include "argparser.h"



int main(int argc, char *argv[]){

    // Process command line arguments.
    auto [parse_result, options] = parseCommandLine(argc, argv); // Parse command line.
    if (parse_result.count("help")) {std:: cout << options.help() << std:: endl; return 0;}
    ParsedValues vals = processParsedValues(parse_result);
    
    // Set up measurement object.
    bool method_type = vals.sampler=="zbaoabz" ? 1 : 0;
    Measurement results = Measurement(method_type, 
                                      vals.burnin, 
                                      vals.t_meas, 
                                      vals.n_dist, 
                                      vals.time_average, 
                                      vals.N_iteration, 
                                      vals.output_name);
    

    // Set up simulation object.
    Simulation simu(vals.sampler, 
                    vals.stepsize, 
                    vals.temperature, 
                    vals.friction, 
                    vals.alpha1, 
                    vals.alpha2, 
                    vals.N_iteration, 
                    vals.t_meas, 
                    vals.burnin, 
                    vals.forcefield, 
                   {vals.init_position[0],vals.init_position[1]}, 
                   {vals.init_velocity[0],vals.init_velocity[1]},
                    results);
    
    // Run simulation.
    simu.run_MPI_simulation(argc, argv);


    return 0;

}



    





