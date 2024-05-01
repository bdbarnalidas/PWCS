OBJS = Vector.o Simulate.o Particle.o Protein.o Force.o Cell.o General.o Collision.o Boundary.o Verlet.o Main.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG) 
LFLAGS = -Wall $(DEBUG)

WholeCellSim : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o WholeCellSim	
	
Vector.o : Vector.cpp Vector.h
	$(CC) $(CFLAGS) Vector.cpp

Simulate.o : Simulate.cpp Simulate.h
	$(CC) $(CFLAGS) Simulate.cpp

Particle.o : Particle.cpp Particle.h
	$(CC) $(CFLAGS) Particle.cpp

Protein.o : Protein.cpp Protein.h Particle.h
	$(CC) $(CFLAGS) Protein.cpp

Force.o : Force.cpp Force.h
	$(CC) $(CFLAGS) Force.cpp

Cell.o : Cell.cpp Cell.h String.h Force.h
	$(CC) $(CFLAGS) Cell.cpp

General.o : General.cpp General.h Protein.h
	$(CC) $(CFLAGS) General.cpp

Collision.o : Collision.cpp Collision.h Vector.h
	$(CC) $(CFLAGS) Collision.cpp

Boundary.o : Boundary.cpp Boundary.h
	$(CC) $(CFLAGS) Boundary.cpp

Verlet.o : Verlet.cpp Verlet.h Force.h
	$(CC) $(CFLAGS) Verlet.cpp
	
Main.o : Main.cpp Simulate.h Cell.h General.h Collision.h Boundary.h Verlet.h
	$(CC) $(CFLAGS) Main.cpp

clean: 
	\rm *.o *~ WholeCellSim
