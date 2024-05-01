#ifndef GENERAL_H
#define GENERAL_H

#include <vector>

#include "Cell.h"

void write_to_xyz(std::ofstream&,Cell&); // write to xyz file
void write_to_xyz(std::ofstream&,std::vector<Cell>&,const Cell&); // write to xyz file

#endif
