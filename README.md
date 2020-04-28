# Solar Sail Controls

A PD controller for a solar sail with controller gains tuned through simulated annealing.

## Getting Started

Run the Main.m file.  This file runs the simulated annealing function on the simulation of a single orbit by calling the simulation.m file in each iteration.  The simulation.m file takes PD controller gains as inputs and calculates the changes in the orbital parameters due to solar radiation pressure.  The effect of the solar radiation pressure varies with the controller gains.  The goal of the simulated annealing algorithm is to determine the optimal controller gains for an orbit raising maneuver.

### Prerequisites

In order to run this simulation, you must have MATLAB installed.


## Built With

* [MATLAB](https://www.mathworks.com/products/matlab.html)


## Authors

* **Kayla Ormiston** 
