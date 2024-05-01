#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Simulate.h"
#include "Cell.h"
#include "General.h"
#include "Collision.h"
#include "Boundary.h"
#include "Verlet.h"

int main()
{
 std::ofstream xyz;
 xyz.open("Trajectory1.xyz");
 Simulate simulate;
 Collision collision(0);
 Cell cell;
 cell.addprotein(cell,1000,2,simulate);
 Boundary boundary(cell.getbreadth(),cell.getlength(),cell.getwidth());
 write_to_xyz(xyz,cell);

 /*for(int i=0;i<simulate.get_nsteps();i++)
 {
  Verlet().VelocityVerlet(cell,cell.getproteins(),simulate,boundary,collision);
  write_to_xyz(xyz,cell);
 }*/
 xyz.close();
}
