# Tetrahedron Packing

### Running the code

Required Installations/Dependencies:  
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
  -b, --bias            induce nearest neighbor attraction bias between tetrahedra to
                        avoid dogpile clustering  
  -np, --noplot         disable plotting tetrahedra with plotly  
  
The default parameters are the parameters I found to produce the best convergence behavior while testing my code.  

### File descriptions
main.py - This file takes in user specified paramaters and runs the packing algorithm with them. After running the algorithm, it prints the basis vectors and vertex coordinates, saves them to a text file, then saves a plotly graph to an html file.  
  
util.py - This file contains the implementation of the packing algorithm itself along with the implementations for a few utility functions, such as a tetrahedron plotting function.  
  
convexPolygon.py - This file contains class definitions for d-dimensional polygons. Because this project deals specifically with tetrahedron packing the file only contains one class, Simplex. As more classes are added, I would consider adding an AbstractPolygon class to enforce consistency and functionality requirements.  
  
tetrahedron-intersect.clj - This file contains postspectacular's clojure implementation of tetrahedron intersection detection, which was linked to in the original problem statement email. To avoid re-implementing intersection detection, I simply integrated this code into my project.  

### Approach
My initial idea to solve this tetrahedron packing problem was to take a gradient descent based approach where n tetrahedra are initialized with random starting positions in space, and then, using a loss function such as container volume or net distance between tetrahedra, their positions are shifted with gradient updates until a minimum is reached. However, some issues I encountered with this approach include the high likelihood to get caught in local minima, the high likelihood to get caught between the gradient and collision constraints, and difficulty formalizing the gradient update rule such that it promotes rotations into the ideal orientations. Inspired by simulated annealing, my final solution attempts to circumvent these issues by updating tetrahedra positions with random translations and rotations that are loosely guided by the goal of minimizing the distance between each tetrahedron and a "center of mass" that is calculated as the elementwise average of all tetrahedron centroid coordinates. This approach is similar in essence to the Monte Carlo simulation described in Chen et al.'s *Dense Crystalline Dimer Packings of Regular Tetrahedra*, except that my move probabilities are not based on an energy calculation, and instead of compressing the tetrahedra with container rescalings, my compression occurs as a result of the distance minimization objective.  
  
Within util.py there are three implementations of the randomizedGuidedPacking described above. The original and V2 implementations are almost identical, except V2 uses a less restrictive rotation update. V3 uses V2 as a base and adds a nearest neighbor translation bias that suggests a more "neat" packing cluster (without the bias, packing many tetrahedra usually results in a cluster with many gaps as the tetrahedra converge to a single point and collide with each other).  

### Conclusions
In general, my solution does not converge to known optimal packings, but it can usually produce a "decent" packing within a few trials. After testing my algorithm with a variety of different parameters, I believe the primary weakness of my approach lies within its simplistic convergence objective. Having the tetrahedra converge to a single point can very quickly produce a rough packing with gaps, but with only this objective, the algorithm struggles to encode the complex orientation requirements of dense packings. I found that the addition of a nearest neighbor bias can reduce the number of, or at least size of, gaps in the final packing, but it still does not produce the densest possible configurations, as dense organizations are not always comprised of tetrahedra oriented symmetrically about the "center of mass". Future work on this project might involve adding a more sophisticated convergence objective to the approach described above.  

