#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <string>
#include <iostream>
#include <cxxopts.hpp>

// DEFAULT PARAMETERS
const std:: string  _N_particles_default {"1000"};
const std:: string  _boxlength_default {"10"};
const std:: string  _forcefield_default {"gauss"};
const std:: string  _init_mode_default {"uniform"};
const std:: string  _beta_default {"10"};
const std:: string  _gamma_default {"0.1"};
const std:: string  _N_iter_default{"10000"};
const std:: string  _stepsize_default {"0.01"};
const std:: string  _N_meas_default {"10"};
const std:: string  _integrator_default {"BAOAB"};
const std:: string  _threads_default {"4"};
const std:: string  _seed_default {"1"};
const std:: string  _output_name_default {"results.csv"};


std:: pair <cxxopts:: ParseResult, cxxopts:: Options> parseCommandLine(int argc, char* argv[]) {
std::string description = R"(
    To run a simulation use IPS.exe and use the flags below to control the settings.
    The program will print a .csv file containing time series data of the mean center of mass
    distance, the mean-squared-displacement, and the kinetic temperature.

    If the --traj flag is passed, the whole trajectory (i.e. all positions of all particles in time)
    is printed to a separate file traj.csv.
    )";
    cxxopts::Options options("IPS.exe", description);

    // Define command line options.
    options.add_options()
        ("N_particles", "Number of particles.",                                         cxxopts:: value <int>()->default_value(_N_particles_default))
        ("boxlength",   "Length of edge of square simulation box.",                     cxxopts:: value <double>()->default_value(_boxlength_default))
        ("forcefield",  "Forcefield between two particles. Either 'gauss' or 'morse'.", cxxopts:: value <std:: string>()->default_value(_forcefield_default))
        ("init_mode",   "Initial positions of the system. Either 'uniform' or 'grid'.", cxxopts:: value <std:: string>()->default_value(_init_mode_default))
        ("beta",        "Inverse temperature parameter in Langevin dynamics.",          cxxopts:: value <double>()->default_value(_beta_default))
        ("gamma",       "Friction parameter in Langevin dynamics.",                     cxxopts:: value <double>()->default_value(_gamma_default))
        ("N_iter",      "Number of simulation steps.",                                  cxxopts:: value <int>()->default_value(_N_iter_default))
        ("stepsize",    "Simulation stepsize.",                                         cxxopts:: value <double>()->default_value(_stepsize_default))
        ("N_meas",      "Take a measurement any 'N_meas' iterations.",                  cxxopts:: value <int>()->default_value(_N_meas_default))
        ("integrator",  "Integrator to be used. Either 'BAOAB' or 'UBU'.",              cxxopts:: value <std:: string>()->default_value(_integrator_default))
        ("threads",     "Number of threads used in the force evaluation.",              cxxopts:: value <int>()->default_value(_threads_default))
        ("seed",        "Randomseed.",                                                  cxxopts:: value <int>()->default_value(_seed_default))
        ("output_name", "Name of the printed file holding the results.",                cxxopts:: value <std:: string>()->default_value(_output_name_default))               
        ("trajectory",  "If flag is set, trajetory will be printed to file.",           cxxopts:: value<bool>())           
        ("help",        "Print help");

    // Parse command line.
    return {options.parse(argc, argv), options};
}


struct ParsedValues{
    int N_particles;
    double boxlength;
    std:: string forcefield;
    std:: string init_mode;
    double beta;
    double gamma;
    int N_iter;
    double stepsize;
    int N_meas;
    std:: string integrator;
    int threads;
    int seed;
    std:: string output_name;
    bool trajectory;
};


ParsedValues processParsedValues(const cxxopts:: ParseResult& result) {
    
    // Access parsed values.
    ParsedValues values;

    values.N_particles = result.count("N_particles") ?   result["N_particles"].as<int>()        : std:: stoi(_N_particles_default);
    values.boxlength   = result.count("boxlength")   ?   result["boxlength"].as<double>()       : std:: stod(_boxlength_default);
    values.forcefield  = result.count("forcefield")  ?   result["forcefield"].as<std::string>() : _forcefield_default;
    values.init_mode   = result.count("init_mode")   ?   result["init_mode"].as<std::string>()  : _init_mode_default;
    values.beta        = result.count("beta")        ?   result["beta"].as<double>()            : std:: stod(_beta_default);
    values.gamma       = result.count("gamma")       ?   result["gamma"].as<double>()           : std:: stod(_gamma_default);
    values.N_iter      = result.count("N_iter")      ?   result["N_iter"].as<int>()             : std:: stoi(_N_iter_default);
    values.stepsize    = result.count("stepsize")    ?   result["stepsize"].as<double>()        : std:: stod (_stepsize_default);
    values.N_meas      = result.count("N_meas")      ?   result["N_meas"].as<int>()             : std:: stoi(_N_meas_default);
    values.integrator  = result.count("integrator")  ?   result["integrator"].as<std::string>() : _integrator_default;
    values.threads     = result.count("threads")     ?   result["threads"].as<int>()            : std:: stoi(_threads_default);
    values.seed        = result.count("seed")        ?   result["seed"].as<int>()               : std:: stoi(_seed_default);
    values.output_name = result.count("output_name") ?   result["output_name"].as<std::string>(): _output_name_default;
    values.trajectory  = result.count("trajectory");

    return values;
}


#endif // ARGPARSER_H