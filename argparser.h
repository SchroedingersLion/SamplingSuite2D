#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <string>
#include <iostream>
#include <cxxopts.hpp>

// DEFAULT PARAMETERS
const std:: string  _sampler_default {"baoab"};
const std:: string  _stepsize_default {"0.1"};
const std:: string  _temperature_default {"1"};
const std:: string  _friction_default {"0.1"};
const std:: string  _alpha1_default {"1"};
const std:: string  _alpha2_default {"1"};
const std:: string  _forcefield_default {"mullerbrown"};
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
    
    cxxopts::Options options("SamplingSuited2D", description);

    // Define command line options.
    options.add_options()
        ("sampler",         "Sampler run ('baoab' or 'zbaoabz').",                      cxxopts:: value <std:: string>()->default_value(_sampler_default))
        ("stepsize",        "Simulation stepsize.",                                     cxxopts:: value <double>()->default_value(_stepsize_default))
        ("temperature",     "Temperature parameter in Langevin dynamics.",              cxxopts:: value <double>()->default_value(_temperature_default))
        ("friction",        "Friction parameter in Langevin dynamics.",                 cxxopts:: value <double>()->default_value(_friction_default))
        ("alpha1",          "Hyperparameter alpha_1 (only used by ZBAOABZ)",            cxxopts:: value <double>()->default_value(_alpha1_default))
        ("alpha2",          "Hyperparameter alpha_2 (only used by ZBAOABZ)",            cxxopts:: value <double>()->default_value(_alpha2_default))        
        ("forcefield",      "Forcefield. Pick 'mullerbrown' or 'morse'.",               cxxopts:: value <std:: string>()->default_value(_forcefield_default))         
        ("init_position",   "Initial position. Enter two comma-separated floats.",      cxxopts:: value<std::vector<double>>()->default_value(_init_position_default))
        ("init_velocity",   "Initial velocity. Enter two comma-separated floats.",      cxxopts:: value<std::vector<double>>()->default_value(_init_velocity_default))
        ("N_iteration",     "Number of simulation steps.",                              cxxopts:: value <int>()->default_value(_N_iteration_default))
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
    double alpha1;
    double alpha2;
    std:: string forcefield;
    std:: vector <double> init_position;
    std:: vector <double> init_velocity;
    int N_iteration;
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
    values.alpha1        = result.count("alpha1")        ? result["alpha1"].as<double>()                     : std:: stod(_alpha1_default);
    values.alpha2        = result.count("alpha2")        ? result["alpha2"].as<double>()                     : std:: stod(_alpha2_default);
    values.forcefield    = result.count("forcefield")    ? result["forcefield"].as<std:: string>()           : _forcefield_default;
    values.init_position = result.count("init_position") ? result["init_position"].as<std::vector<double>>() : std:: vector<double> {0,0};
    values.init_velocity = result.count("init_velocity") ? result["init_velocity"].as<std::vector<double>>() : std:: vector<double> {0,0};
    values.N_iteration   = result.count("N_iteration")   ? result["N_iteration"].as<int>()                   : std:: stoi(_N_iteration_default);
    values.t_meas        = result.count("t_meas")        ? result["t_meas"].as<int>()                        : std:: stoi(_t_meas_default);
    values.n_dist        = result.count("n_dist")        ? result["n_dist"].as<int>()                        : std:: stoi(_n_dist_default);
    values.burnin        = result.count("burnin")        ? result["burnin"].as<int>()                        : std:: stoi(_burnin_default);
    values.seed          = result.count("seed")          ? result["seed"].as<int>()                          : std:: stoi(_seed_default);
    values.output_name   = result.count("output_name")   ? result["output_name"].as<std::string>()           : _output_name_default;
    values.time_average  = result.count("time_average");

    return values;
}


#endif // ARGPARSER_H