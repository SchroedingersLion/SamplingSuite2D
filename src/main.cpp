

#include "parameters.h"
#include "measurement.h"
#include "simulation.h"   
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
                    vals.Sundman_m,
                    vals.Sundman_M, 
                    vals.N_iteration, 
                    vals.t_meas, 
                    vals.burnin,
                    vals.seed, 
                    vals.forcefield, 
                   {vals.init_position[0],vals.init_position[1]}, 
                   {vals.init_velocity[0],vals.init_velocity[1]},
                    vals.sigma_noise,
                    results);
    
    // Run simulation.
    simu.run_MPI_simulation(argc, argv);


    return 0;

}



    





