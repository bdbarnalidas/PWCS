#ifndef CELL_H
#define CELL_H

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "Protein.h"
#include "Simulate.h"
//#include "Boundary.h"

class Cell
{
 protected:
 double length; // measurements along the Y-axis
 double breadth; // measurements along the X-axis
 double width; // measurements along the Z-axis
 int n_protein; // total no of proteins present inside the cell
 std::vector <Protein> proteins; // vector proteins will store all the proteins in the cell
 //Boundary boundary; // the cell's boundary
 double row_no; // row no of the subcell in the entire cell
 double col_no; // col no of the subcell in the entire cell
 double startx; // starting point along X-axis
 double starty; // starting point along Y-axis
 double startz; //starting point along Z-axis
 double endx; // ending point along X-axis
 double endy; // ending point along Y-axis
 double endz; // ending point along Z-axis
 public:
 Cell(); // default constructor
 Cell(const double,const double,const double,const int,const double,const double,const double,const double); // overloading constructor
 double getlength() const; // to get the value of length
 double getbreadth() const; // to get the value of breadth
 double getwidth() const; // to get the value of width
 int getprotein() const; // to get the no of proteins inside the cell
 void print() const; // to print all the parameters of the cell
 void addprotein(const Cell&,const int,const int,const Simulate&); // to add proteins in the cell (modify member variable -> proteins)
 std::vector<Protein>& getproteins(); // return the proteins vector of cell
 //Boundary& get_boundary(); // return the boundary of the cell
 double get_rowno() const; // return the row no of the subcell
 double get_colno() const; // return the col no of the subcell
 void set_rowno(const double); // set the row_no of the subcell
 void set_colno(const double); // set the col_no of the subcell
 double get_startx() const; // return the value of startx
 double get_starty() const; // return the value of starty
 double get_startz() const; // return the value of startz
 void set_startx(const double); // set the value of startx
 void set_starty(const double); // set the value of starty
 void set_startz(const double); // set the value of startz
 double get_endx() const; // return the value of endx
 double get_endy() const; // return the value of endy
 double get_endz() const; // return the value of endz
 void set_endx(const double); // set the value of endx
 void set_endy(const double); // set the value of endy
 void set_endz(const double); // set the value of endz
};

#endif
