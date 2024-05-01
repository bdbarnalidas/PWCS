#include "General.h"

#include "Protein.h"

void write_to_xyz(std::ofstream& xyz,std::vector<Cell>& subcells) // write to xyz file
{
 int size=subcells.size();
 int count_protein=0;
 for(int i=0;i<size;i++)
 count_protein+=subcells[i].getprotein(); // find out the total no of proteins in the entire cell
 xyz<<count_protein<<"\n";
 for(int i=0;i<size;i++)
 {
  std::vector<Protein>& pro=subcells[i].getproteins();
  int number=pro.size();
  for(int j=0;j<number;j++)
  xyz<<pro[j].getsymbol()<<" "<<(pro[j].getposition()).getx()<<" "<<(pro[j].getposition()).gety()<<" "<<(pro[j].getposition()).getz()<<"\n";
 }
}

void write_to_xyz(std::ofstream& xyz,Cell& cell) // write to xyz file
{
 xyz<<cell.getprotein()<<"\n";
 std::vector<Protein>& pro=cell.getproteins();
 int number=pro.size();
 for(int j=0;j<number;j++)
 xyz<<pro[j].getsymbol()<<" "<<(pro[j].getposition()).getx()<<" "<<(pro[j].getposition()).gety()<<" "<<(pro[j].getposition()).getz()<<"\n";
}
