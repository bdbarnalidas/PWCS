#ifndef SIMULATE_H
#define SIMULATE_H

#include <string>
#include <sstream>
#include <fstream>

class Simulate
{
 private:
 double duration; // total duration of the simulation
 double timeInterval; // delta_t
 int nsteps; // no. of steps for running the simulation
 double type; // type of force-field = attractive/repulsive
 int nrow; // total no of rows in the cell
 int ncol; // total no of cols in the cell
 public:
 Simulate(); // default constructor
 Simulate(const double,const int,const int,const int); // overloading constructor
 double get_timeInterval() const; // read the value of timeInterval
 int get_nsteps() const; // get the number of steps
 double get_type() const; // get the type of force-field
 int get_nrow() const; // get the no of rows
 int get_ncol() const; // get the no of cols
};

#endif
