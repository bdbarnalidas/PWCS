# PWCS

High-performance whole-cell simulation exploiting modular cell biology principles

We provide CUDA C/C++ simulation code for performing parallel whole-cell simulation exploiting the existing HPC architecture. This includes an efficient simulation guided by our hashing based data structure called the Cellular Dictionary. Important features of this simulation tool are summarized below.

    Spatially localized densely connected protein clusters/modules directed simulation.
    Can properly handle boundary conditions.
    Can detect and resolve collision scenarios.


Please download and unzip 'Executable.zip'.

The extracted folder 'Executable' consists of 2 linux executable files 'Serial_Simulator' and 'Parallel_Simulator'.  

******************************
Serial_Simulator (CPU_Version)
******************************

Serial_Simulator can perform whole-cell simulation of the unicellular bacterium, Escherichia coli, on CPU (1 core).

Use the command './Serial_Simulator' for executing the serial simulation.

********************************
Parallel_Simulator (GPU_Version)
********************************

Parallel_Simulator can perform parallel whole-cell simulation of Escherichia coli on 128 GPU cores.

Use the command './Parallel_Simulator' for executing the parallel simulation.

For any doubts or suggestions please contact
barnali.das@iitkgp.ac.in

Please visit https://cosmos.iitkgp.ac.in/PWCS/




