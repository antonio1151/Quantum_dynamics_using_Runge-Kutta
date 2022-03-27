# Quantum_dynamics_using_Runge-Kutta

This is a Fortran program for performing Quantum dynamics for any system. The dynamics is performed by integrating the Time-Dependent Schrodinger equation using Runge-Kutta of order 4th. 

To run it, you need to introduce your initial state in a file named "INITIAL_STATE.dat". In the first column, you must put the values of the real part of the coefficients of the initial state and in the second column the imaginary part. 

To adapt it to your own problem, you must change the Hamiltonian in the subroutine "HH", within the main file "dynamics_rk", for the one you want to study. Additionally, you can play with the timestep, time for printing information, and the final time 


The program comes with an example of two-level system interacting with a laser pulse, where the initial state is the lower energy state "a". The dynamics of the populations of the energy states "a" and "b" are shown in the figure "results.png"


To compile it you have to use the command:

XXXX  dynamics_rk RK4.f

where XXXX is the FORTRAN compiler you use, RK4.f is the subroutine that perform the integration using Runge-Kutta of order 4th.


