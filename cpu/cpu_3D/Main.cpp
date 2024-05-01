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
 xyz.open("Trajectory.xyz");
 Simulate simulate;
 Collision collision(1);
 Cell wholecell;
 srand(time(NULL));
 int limit=simulate.get_nrow()*simulate.get_ncol();
 std::vector <Cell> Cells(limit);
 int i=0; int count=0;
 while(i<simulate.get_nrow())
 {
  int colcount=0;
  while(colcount<simulate.get_ncol())
  {
   Cell cell(wholecell.getlength(),wholecell.getbreadth(),wholecell.getwidth(),(rand()%7+1),i,colcount,(wholecell.getbreadth()+1)*colcount,(wholecell.getlength()+1)*i);
   cell.addprotein(cell,cell.getprotein(),2,simulate);
   Cells[count]=cell;
   colcount++;count++;
  }
  i++;
 }

 for(int i=0;i<limit;i++)
 {
  std::cout<<"Row no= "<<Cells[i].get_rowno()<<" Col no= "<<Cells[i].get_colno()<<"\n";
  std::cout<<"startx= "<<Cells[i].get_startx()<<" starty= "<<Cells[i].get_starty()<<" startz= "<<Cells[i].get_startz()<<"\n";
  std::cout<<"endx= "<<Cells[i].get_endx()<<" endy= "<<Cells[i].get_endy()<<" endz= "<<Cells[i].get_endz()<<"\n";
  std::cout<<"Boundary a= "<<((Cells[i].get_boundary()).geta()).getx()<<" "<<((Cells[i].get_boundary()).geta()).gety()<<" "<<((Cells[i].get_boundary()).geta()).getz()<<"\n";
  std::cout<<"Boundary b= "<<((Cells[i].get_boundary()).getb()).getx()<<" "<<((Cells[i].get_boundary()).getb()).gety()<<" "<<((Cells[i].get_boundary()).getb()).getz()<<"\n";
  std::cout<<"Boundary c= "<<((Cells[i].get_boundary()).getc()).getx()<<" "<<((Cells[i].get_boundary()).getc()).gety()<<" "<<((Cells[i].get_boundary()).getc()).getz()<<"\n";
  std::cout<<"\n\n\n";
 }




 int size=Cells.size();
 for(int j=0;j<size;j++)
 {
  for(int i=0;i<simulate.get_nsteps();i++)
  {
   Verlet().VelocityVerlet(Cells[j],Cells[j].getproteins(),simulate,Cells[j].get_boundary(),collision);
   write_to_xyz(xyz,Cells,wholecell);
  }
 }


 xyz.close();
}
