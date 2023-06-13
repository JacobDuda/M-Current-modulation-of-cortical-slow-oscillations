# M-Current-modulation-of-cortical-slow-oscillations

This repository is part of the paper: M-current modulation of cortical slow oscillations: network dynamics and computational modelling, by: Leonardo Dalla Porta, Almudena Barbero-Castillo, Jose Manuel Sanchez-Sanchez and Maria V. Sanchez-Vives.

Plos Computational Biology: In Press

Whenever this code is used, please cite the aforomentioned paper.

########

The code is implemeted in c++ language.

The code can be compiled as: g++ main.cpp -O3 -march=native -o exec

In order to execute, run: ./exec 1.0 123456
The first parameter is the M-current concentration (a factor multiplying M-current maximal conductance g_M)
123456 is the random number seed.

Tested in: gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0

Contact: leonardodallaporta@gmail.com
