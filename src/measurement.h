#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <vector>
#include <string>
#include <fstream>
#include <mpi.h>
#include "parameters.h"


/* THIS DESCRIPTION NEEDS TO BE UPDATED  
This is the measurement class that is used by the samplers to obtain observables. 
    In order to modify what observables are collected, the user has to do 2 things:
    a) in the constructor, they have to specify the numbers of observables to be taken by
       adjusting the value of the variable "no_observables". 
    b) in the function "take_measurement" the user has to adjust the formulas used to compute an observable from 
       the parameters. */

class Measurement {

    public:

        // constructor
        Measurement (const int method_type, 
                     const int burnin, 
                     const int t_meas, 
                     const int n_dist,
                     const bool time_average, 
                     const int N_iteration,
                     const std:: string output_name):
                     method_type {method_type}, 
                     burnin {burnin}, 
                     t_meas {t_meas}, 
                     n_dist {n_dist},
                     time_average {time_average}, 
                     N_iteration {N_iteration},
                     output_name {output_name}
            {
                
                /*######## ENTER THE NUMBER OF OBSERVABLES TO COLLECT ############*/
                no_observables = 4; 
                /*################################################################*/
                

                observables.resize(no_observables);
                observable_sums.resize(no_observables);
                col_names.resize(no_observables);

                int no_elements {(N_iteration-burnin) / (n_dist*t_meas)}; // Number of elements needed in the t-averaged results vector.
                if (no_elements<0) throw std:: invalid_argument( "\n Combination of N_iteration, burnin, n_dist, and t_meas makes no sense. Is N_iter < burnin? \n" );

                observable_tavgs = std:: vector <std:: vector <float>> (no_observables, std::vector <float> (no_elements));

                // the adaptive schemes also collect the adaptive step sizes.
                if (method_type==1){ 
                    dt_vals_raw.resize(no_elements);
                    zeta_raw.resize(no_elements);
                } 
            
            };



        void take_measurement(const parameters& parameters){

            /* ########### COMPUTE CURRENT OBSERVABLE VALUES FROM PARAMETERS #####################################
               The number of entries in vector "observables" must correspond to member variable "no_observables" set by the user
               in the constructor above. The reweighting for the adaptive schemes is done automatically by the "process_sample" routine. 
               The adaptive schemes will also automatically store the adaptive step size and the \zeta variable. */            
            observables[0] = parameters.position.x;	// x-coordinate
            observables[1] = parameters.position.y;  // y-coordinate
            observables[2] = 0.5*(parameters.velocity.x*parameters.velocity.x + parameters.velocity.y*parameters.velocity.y);   // Tkin
            observables[3] = -0.5*(parameters.position.x*parameters.force.x + parameters.position.y*parameters.force.y);        // Tconf
            /*########################################################################*/

            /*################# ENTER NAMES OF OBSERVABLES (WILL BE HEADER OF OUTPUTFILE)####*/
            col_names[0] = "x";
            col_names[1] = "y";
            col_names[2] = "Tkin";
            col_names[3] = "Tconf";
            /*################################################################################*/

            process_sample(parameters);   // compute (reweighted) time average and stores results for later print-out.

            return;

        };


        void mpi_reduction(MPI_Comm& comm, const int& rank, const int& nr_proc);  // for mpi average, implemented in .cpp file.
        void print_to_csv();                         // for print out, implemented in .cpp file.


    private:
        std:: vector <double> observable_sums;          // sum of observable samples for time average
        std:: vector <double> observables;              // vector storing the new sample (one entry per observable)
        std:: vector <std:: string> col_names;          // names of the columns in the output file (names of the observables).
        std:: vector <std:: vector <float>> observable_tavgs;       // store the evolving time average for each observable
        std:: vector <std:: vector <float>> observable_printout;    // store the mpi process average of vector "observable_tavgs" on rank 0 
        std:: vector <double> dt_vals_raw;              // store samples of adaptive step size (only for adaptive schemes)
        std:: vector <double> zeta_raw;                 // store samples of zeta (only for adaptive schemes)                 
        int no_observables;                             // number of observables to be taken
        const int method_type;                          // 0=constant step size scheme, 1=adaptive scheme
        const bool time_average;                        // decide whether observables will be time-averaged
        int k{0}, ctr{0}, t_avg_normalizer{0};          // help variables
        const int burnin, t_meas, n_dist, N_iteration;     
        double dt_sum{0}; 
        const std:: string output_name;                              
        void process_sample(const parameters& parameters);   // add new sample to time average and store result, implemented in .cpp file.

};



// implementation of measurement class routines that the user does not need to modify.

inline void Measurement:: print_to_csv(){
    
    std:: ofstream file{output_name};
    std:: cout << "Writing to file...\n";

    // Write header with specified column names.
    file << "Iteration ";
    for (size_t k=0; k<col_names.size(); ++k){
        file << col_names[k] << " ";
    }
    if (method_type==1) file << "dt zeta "; 
    file << "\n";

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
        
        if (method_type==1) file << dt_vals_raw[i] << zeta_raw[i] << " ";
        
        file << "\n";
    }

    file.close();

    return;

}


inline void Measurement:: mpi_reduction(MPI_Comm& comm, const int& rank, const int& nr_proc){   // MPI AVERAGE ROUTINE
    
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


inline void Measurement:: process_sample(const parameters& parameters){

    if (method_type==0){   // Constant step size scheme, eg. BAOAB.

        // Take care of time average.
        if (time_average){
            for (int i=0; i<no_observables; ++i) observable_sums[i] += observables[i];  // Add to sum for time average.

            ++t_avg_normalizer;

            // Store current moving average value for print-out, but only every n_dist steps.        
            if (ctr % n_dist == 0){                                             
                for (int i=0; i<no_observables; ++i) observable_tavgs[i][k] = observable_sums[i] / t_avg_normalizer;
                ++k;
            }
        }
        // No time average.
        else{
             // Store the observables for print-out without averaging.
             if (ctr % n_dist == 0){                                             
                for (int i=0; i<no_observables; ++i) observable_tavgs[i][k] = observables[i];
                ++k;
             }           
        }          
    }

    else if (method_type==1){   // Adaptive step size scheme.

        if (time_average){
            for (int i=0; i<no_observables; ++i) observable_sums[i] += observables[i] * parameters.dt;  // Reweighting.
            dt_sum += parameters.dt;

            if (ctr % n_dist == 0){                                                 
                for (int i=0; i<no_observables; ++i) observable_tavgs[i][k] = observable_sums[i] / dt_sum;
                dt_vals_raw[k] = parameters.dt;
                zeta_raw[k] = parameters.zeta;
                ++k;
            }
        }
        else {
            if (ctr % n_dist == 0){                                                 
                for (int i=0; i<no_observables; ++i) observable_tavgs[i][k] = observables[i];
                dt_vals_raw[k] = parameters.dt;
                zeta_raw[k] = parameters.zeta;
                ++k;
            }           
        }         
    }

    ++ctr;

}


#endif // MEASUREMENT_H