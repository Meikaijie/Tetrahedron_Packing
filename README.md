# Tetrahedron Packing

### Running the code

Required Installations:  
Python 2,
clojure,
numpy,
plotly  
  
To run with default parameters, simply enter python main.py  
For a description of all parameters, enter python main.py -h  

For convenience, the latter command outputs the following:  
optional arguments:  
  -h, --help            show this help message and exit  
  -n NUMTETRAS, --numtetras NUMTETRAS
                        number of tetrahedra to pack, valid at or above 1  
  -is INITSTEPSCALE, --initstepscale INITSTEPSCALE
                        initial step size coefficient, valid between 0 and 1  
  -sr STEPSCALEREDUCTION, --stepscalereduction STEPSCALEREDUCTION
                        step scale decay rate, valid between 0 and 1  
  -st STEPTHRESHOLD, --stepthreshold STEPTHRESHOLD
                        minimum allowed step scale, valid above 0  
  -ir INITROTATIONRANGE, --initrotationrange INITROTATIONRANGE
                        initial permitted range of a single rotation move  
  -rr ROTREDUCTION, --rotreduction ROTREDUCTION
                        rotation range decay rate, valid between 0 and 1  
  -rt ROTATIONTHRESHOLD, --rotationthreshold ROTATIONTHRESHOLD
                        minimum allowed rotation range  
  -l LENGTH, --length LENGTH
                        edge length of the tetrahedra  
  -i ITERATIONS, --iterations ITERATIONS
                        max number of iterations  
  -v, --verbose  
  -b, --bias            induce pairwise attraction bias between tetrahedra to
                        avoid dogpile clustering  
  -np, --noplot         disable plotting tetrahedra with plotly  
  
### File descriptions
  
### Approach