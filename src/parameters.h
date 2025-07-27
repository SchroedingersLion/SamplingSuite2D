#ifndef PARAMETERS_H
#define PARAMETERS_H

struct coordinate{
    double x{0}, y{0};
};

struct parameters{
    coordinate position, velocity, force;
    double zeta, dt;    // Used only by ZBAOABZ.
};


#endif //PARAMETERS