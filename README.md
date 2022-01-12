# Symplectic Double Pendulum Simulator
Double pendulum simulator using a symplectic Euler's method and Hamiltonian mechanics. The program calculates the momentum and position of each ball using Hamiltonian formulations, and uses those calculations to adjust the pendulum's course. This provides a double pendulum that stays stable and more analagous to real-life for longer periods of time, without the need of more complex math/animation modules. 

<p align="center">
<img width="400" alt="Screen Shot 2022-01-11 at 4 01 07 PM" src="https://user-images.githubusercontent.com/50428986/149040320-1d5251a4-438f-45db-9f3f-0a44928e3f00.png">
<p>

# Functionality
- Visualization UI that shows scaled changes in pendulum lengths, masses, and starting angles
- Generates EPS plot of path and phase diagram on command
- Built to run on a Rasberry Pi/Linux based system
  
# Usage Instructions
Run `./doublependulum.py` in Linux terminal.
If pendulum animation becomes unstable, the momenta of the masses are most likely too large for the allowed framerate of matplotlib. Try lowering masses or lengthening pendulums.
  
# Required Modules
- numpy
- matplotlib
- threading
- time

# Contributors
Scott Marino / scottmarino@icloud.com
