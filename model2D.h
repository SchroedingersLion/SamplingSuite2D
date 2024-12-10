#ifndef MODEL2D
#define MODEL2D

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>


struct coordinate{    // Used to denote positions, velocities and forces.
    double x{0}, y{0};
};


// ###################### Model CLASS DEFINITION ##################################################################################

class Model2D {

    public: 

        // MEMBERS THAT WILL BE ACCESSED BY MEASUREMENT OR SIMULATION CLASSES
        coordinate position, velocity, force;
        coordinate (Model2D::* get_force) (); // Points to force function.

        // CONSTRUCTOR.
        Model2D(const std:: string& forcefield, 
                 const coordinate init_position, 
                 const coordinate init_velocity)
                : forcefield {forcefield},
                  position {init_position},
                  velocity {init_velocity} 
            {

                std:: cout  << "Building 2D model " << forcefield << std:: endl;

                // Specify force field.
                if (forcefield=="MullerBrown") get_force = &Model2D:: get_force_MullerBrown;
                else throw std:: invalid_argument( "Invalid forcefield. See --help." );

            }

    private:
        const std:: string forcefield;  // Specifies interaction potential.
        

        // ########## FORCES ###############################################################################################################        
        coordinate get_force_MullerBrown();            // Muller-Brown interaction.
        //##################################################################################################################################

};
// ##################### END OF CLASS DEFINITION ##############################################







// ####### DEFINE FORCE FUNCTIONS ######################

// constants used in the force routine
const std:: vector <double> MullerBrown_A {-200, -100, -170, 15};  
const std:: vector <double> MullerBrown_a {-1, -1, -6.5, 0.7};
const std:: vector <double> MullerBrown_b {0, 0, 11, 0.6};
const std:: vector <double> MullerBrown_c {-10, -10, -6.5, 0.7};
const std:: vector <double> MullerBrown_x {1, 0, -0.5, -1};
const std:: vector <double> MullerBrown_y {0, 0.5, 1.5, 1};
double xdiff, ydiff, exponent;


inline coordinate  Model2D:: get_force_MullerBrown(){

    force.x = force.y = 0;

    for (int i=0; i<4; ++i){
        xdiff = position.x - MullerBrown_x[i];
        ydiff = position.y - MullerBrown_y[i];
        exponent = MullerBrown_a[i]*xdiff*xdiff + MullerBrown_b[i]*xdiff*ydiff + MullerBrown_c[i]*ydiff*ydiff;

        force.x -= MullerBrown_A[i] * exp( exponent ) * (2*MullerBrown_a[i]*xdiff + MullerBrown_b[i]*ydiff);
        force.y -= MullerBrown_A[i] * exp( exponent ) * (2*MullerBrown_c[i]*ydiff + MullerBrown_b[i]*xdiff);
    }

};

// ########### END OF MEMBER DEFINITIONS ##############################################


#endif // 2D_MODEL


