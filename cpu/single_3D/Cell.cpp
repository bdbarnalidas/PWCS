#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Cell.h"
#include "String.h"
#include "Force.h"

Cell::Cell() // default constructor
{
 length=200;
 breadth=200;
 width=200;
 n_protein=1000;
 row_no=0;
 col_no=0;
 startx=0;
 starty=0;
 startz=0;
 endx=breadth;
 endy=length;
 endz=width;
 //boundary.seta(length,0,0);
 //boundary.setb(0,breadth,0);
 //boundary.setc(0,0,width);
}

Cell::Cell(const double a,const double b,const double c,const int d,const double e,const double f,const double g,const double h) // overloading constructor
{
 length=a;
 breadth=b;
 width=c;
 n_protein=d;
 row_no=e;
 col_no=f;
 startx=g;
 starty=h;
 startz=0;
 endx=startx+breadth;
 endy=starty+length;
 endz=width;
 //boundary.seta(endx,(length+1)*row_no,startz);
// boundary.setb((breadth+1)*col_no,endy,startz);
 //boundary.setc((breadth+1)*col_no,(length+1)*row_no,endz);
}

double Cell::getlength() const // read the value of length
{
 return this->length;
}

double Cell::getbreadth() const // read the value of breadth
{
 return this->breadth;
}

double Cell::getwidth() const // read the value of width
{
 return this->width;
}

int Cell::getprotein() const // read the total no of proteins inside the cell
{
 return n_protein;
}

void Cell::print() const // print the values of all the parameters of the objects of the Cell class
{
 std::cout<<"\n\nLength="<<length<<" Breadth="<<breadth<<" Width="<<width<<"\n";
 std::cout<<"Total no of proteins in the cell= "<<n_protein<<"\n\n";
}

void Cell::addprotein(const Cell& obj,const int n,const int seed,const Simulate& sim) // add proteins to cell and initializing the parameters of protein
{
 double a,b,c; // variables to store the randomly generated coordinates
 std::string symbol; // string to store the symbol of the protein
 srand(seed); // randomize seed
 for(int i=0;i<n;i++)
 {
  a=((double) (rand() % 99 + 1) / 100.0) * obj.getlength();
  b=((double) (rand() % 99 + 1) / 100.0) * obj.getbreadth();
  c=((double) (rand() % 99 + 1) / 100.0) * obj.getwidth();
  symbol="Prot"+patch::to_string(i); // symbol of protein
  Vector position(a,b,c); // position vector
  Vector velocity(0,0,0); // velocity vector, initially at rest
  Vector halfvelocity(0,0,0); // halfvelocity vector, initially at rest
  Vector force(0,0,0); // force vector, not calculated yet
  if(i==10)
  {
   symbol="Barn"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,2.0,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else if(i==20)
  {  
   symbol="C"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,2.5,position,velocity,halfvelocity,force); // create a protein 
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else if (i==30)
  {
   symbol="H"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,1.5,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else if (i==40)
  {
   symbol="Z"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,1.0,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else if (i==50)
  {
   symbol="Cu"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,1.8,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else if (i==60)
  {
   symbol="O"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,0.9,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else if (i==70)
  {
   symbol="N"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,2.2,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else if (i==80)
  {
   symbol="S"+patch::to_string(i); // symbol of protein
   Protein p(symbol,1,1.0,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
  else
  {
   Protein p(symbol,1,0.7,position,velocity,halfvelocity,force); // create a protein
   proteins.push_back(p); // push the protein into the vector proteins
  }
 }
 Force().GravitationalForces(proteins,sim); // compute the forces acting on each protein by using the force-field
}

std::vector<Protein>& Cell::getproteins() // return the proteins vector of the cell
{
 return proteins;
}

/*Boundary& Cell::get_boundary() // get the boundary of the cell
{
 return boundary;
}*/

double Cell::get_rowno() const // get the row_no of the subcell
{
 return this->row_no;
}

double Cell::get_colno() const // get the col_no of the subcell
{
 return this->col_no;
}

void Cell::set_rowno(const double a) // set the row_no of the subcell
{
 this->row_no=a;
}

void Cell::set_colno(const double a) // set the col_no of the subcell
{
 this->col_no=a;
}

double Cell::get_startx() const // get the value of startx
{
 return this->startx;
}

double Cell::get_starty() const // get the value of starty
{
 return this->starty;
}

double Cell::get_startz() const // get the value of startz
{
 return this->startz;
}

void Cell::set_startx(const double a) // set the value of startx
{
 this->startx=a;
}

void Cell::set_starty(const double a) // set the value of starty
{
 this->starty=a;
}

void Cell::set_startz(const double a) // set the value of startz
{
 this->startz=a;
}

double Cell::get_endx() const // get the value of endx
{
 return this->endx;
}

double Cell::get_endy() const // get the value of endy
{
 return this->endy;
}

double Cell::get_endz() const // get the value of endz
{
 return this->endz;
}

void Cell::set_endx(const double a) // set the value of endx
{
 this->endx=a;
}

void Cell::set_endy(const double a) // set the value of endy
{
 this->endy=a;
}

void Cell::set_endz(const double a) // set the value of endz
{
 this->endz=a;
}

