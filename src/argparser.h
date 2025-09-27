#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <string>
#include <iostream>
#include "cxxopts.hpp"

// DEFAULT PARAMETERS.
const std:: string  _sampler_default {"baoab"};
const std:: string  _stepsize_default {"0.1"};
const std:: string  _temperature_default {"1"};
const std:: string  _friction_default {"0.1"};
const std:: string  _alpha_default {"1"};
const std:: string  _omega_default {"1"};
const std:: string  _Sundman_m_default {"0.1"};
const std:: string  _Sundman_M_default {"10"};
const std:: string  _forcefield_default {"muellerbrown"};
const std:: string  _init_position_default {"0,0"};
const std:: string  _init_velocity_default {"0,0"};
const std:: string  _N_iteration_default{"10000"};
const std:: string  _t_meas_default {"10"};
const std:: string  _n_dist_default {"1"};
const std:: string  _burnin_default {"0"};
const std:: string  _seed_default {"1"};
const std:: string  _output_name_default {"results.csv"};


std:: pair <cxxopts:: ParseResult, cxxopts:: Options> parseCommandLine(int argc, char* argv[]) {
    
    std::string description = R"(
        To run a simulation run "mpirun -n N SamplingSuite2D" and add the flags below to control the settings. N denotes the number of mpi processes.
        
        The program will print a .csv file containing time series data of the observables specified in measurement.h.

        To run multiple trajectories, run the code with multiple mpi processes. The output file will contain the inter-trajectory averages.
        To also time-average the measurements, pass the corresponding flag (see below).

        The first column in the output file will always contain iteration count. The next columns hold the observables. In case of ZBAOABZ, the last two columns
        will hold values of \zeta and the adaptive stepsize dt. These two columns will NEVER be averaged. In case of multiple trajectories, the values in 
        these columns stem from the first trajectory.

        )";
    
    cxxopts::Options options("SamplingSuite2D", description);

    // Define command line options.
    options.add_options()
        ("sampler",         "Sampler to run ('euler', 'baoab' or 'zbaoabz').",          cxxopts:: value <std:: string>()->default_value(_sampler_default))
        ("stepsize",        "Simulation stepsize.",                                     cxxopts:: value <double>()->default_value(_stepsize_default))
        ("temperature",     "Temperature parameter in Langevin dynamics.",              cxxopts:: value <double>()->default_value(_temperature_default))
        ("friction",        "Friction parameter in Langevin dynamics.",                 cxxopts:: value <double>()->default_value(_friction_default))
        ("alpha",           "Hyperparameter alpha (only used by ZBAOABZ)",              cxxopts:: value <double>()->default_value(_alpha_default))
        ("omega",           "Hyperparameter omega (only used by ZBAOABZ)",              cxxopts:: value <double>()->default_value(_omega_default))
        ("sundman_m",       "Hyperparameter m (only used by ZBAOABZ)",                  cxxopts:: value <double>()->default_value(_Sundman_m_default))
        ("sundman_M",       "Hyperparameter M (only used by ZBAOABZ)",                  cxxopts:: value <double>()->default_value(_Sundman_M_default))           
        ("forcefield",      "Forcefield. Allowed values are 'muellerbrown', 'doublewell', 'ackley'," 
                            "'rosenbrock', 'beale', 'entropicbarrier'," 
                            "'star', 'neal', or 'harmonic_asym'.",                      cxxopts:: value <std:: string>()->default_value(_forcefield_default))         
        ("init_position",   "Initial position. Enter two comma-separated floats.",      cxxopts:: value<std::vector<double>>()->default_value(_init_position_default))
        ("init_velocity",   "Initial velocity. Enter two comma-separated floats.",      cxxopts:: value<std::vector<double>>()->default_value(_init_velocity_default))
        ("N_iteration",     "Number of simulation steps.",                              cxxopts:: value <long long>()->default_value(_N_iteration_default))
        ("t_meas",          "Take a measurement any 't_meas' iterations.",              cxxopts:: value <int>()->default_value(_t_meas_default))
        ("n_dist",          "Print any n_dist TAKEN measurements to output file.",      cxxopts:: value <int>()->default_value(_n_dist_default))
        ("burnin",          "Start taking measurements only after burnin iterations.",  cxxopts:: value <int>()->default_value(_burnin_default))
        ("seed",            "Randomseed.",                                              cxxopts:: value <int>()->default_value(_seed_default))
        ("output_name",     "Name of the printed file holding the results.",            cxxopts:: value <std:: string>()->default_value(_output_name_default))               
        ("time_average",    "If flag is set, results will be time-averaged on the fly.",cxxopts:: value<bool>())       
        ("help",            "Print help");

    // Parse command line.
    return {options.parse(argc, argv), options};
}


struct ParsedValues{
    std:: string sampler;
    double stepsize;
    double temperature;
    double friction;
    double alpha;
    double omega;
    double Sundman_m;
    double Sundman_M;
    std:: string forcefield;
    std:: vector <double> init_position;
    std:: vector <double> init_velocity;
    long long N_iteration;
    int t_meas;
    int n_dist;
    int burnin;
    int seed;
    std:: string output_name;
    bool time_average;
};


ParsedValues processParsedValues(const cxxopts:: ParseResult& result) {
    
    // Access parsed values.
    ParsedValues values;

    values.sampler       = result.count("sampler")       ? result["sampler"].as<std:: string>()              : _sampler_default;
    values.stepsize      = result.count("stepsize")      ? result["stepsize"].as<double>()                   : std:: stod(_stepsize_default);
    values.temperature   = result.count("temperature")   ? result["temperature"].as<double>()                : std:: stod(_temperature_default);
    values.friction      = result.count("friction")      ? result["friction"].as<double>()                   : std:: stod(_friction_default);
    values.alpha         = result.count("alpha")         ? result["alpha"].as<double>()                      : std:: stod(_alpha_default);
    values.omega         = result.count("omega")         ? result["omega"].as<double>()                      : std:: stod(_omega_default);
    values.Sundman_m     = result.count("sundman_m")     ? result["sundman_m"].as<double>()                  : std:: stod(_Sundman_m_default); 
    values.Sundman_M     = result.count("sundman_M")     ? result["sundman_M"].as<double>()                  : std:: stod(_Sundman_M_default); 
    values.forcefield    = result.count("forcefield")    ? result["forcefield"].as<std:: string>()           : _forcefield_default;
    values.init_position = result.count("init_position") ? result["init_position"].as<std::vector<double>>() : std:: vector<double> {0,0};
    values.init_velocity = result.count("init_velocity") ? result["init_velocity"].as<std::vector<double>>() : std:: vector<double> {0,0};
    values.N_iteration   = result.count("N_iteration")   ? result["N_iteration"].as<long long>()             : std:: stoi(_N_iteration_default);
    values.t_meas        = result.count("t_meas")        ? result["t_meas"].as<int>()                        : std:: stoi(_t_meas_default);
    values.n_dist        = result.count("n_dist")        ? result["n_dist"].as<int>()                        : std:: stoi(_n_dist_default);
    values.burnin        = result.count("burnin")        ? result["burnin"].as<int>()                        : std:: stoi(_burnin_default);
    values.seed          = result.count("seed")          ? result["seed"].as<int>()                          : std:: stoi(_seed_default);
    values.output_name   = result.count("output_name")   ? result["output_name"].as<std::string>()           : _output_name_default;
    values.time_average  = result.count("time_average");

    return values;
}


#endif // ARGPARSER_H