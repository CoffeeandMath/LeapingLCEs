# LeapingLCEs
Modeling deformations of a LCE sheet

## Installation
This code relies on Deal.II to do finite element calculations. For detailed instructions and documentation on this software, refer to https://www.dealii.org/

After Deal.II is installed, navigate to the source directory (on the level of the license and README files) and enter the following commands:


```ruby
	mkdir build && cd build/
	cmake ..
	make
	mkdir solutions
```


The first three commands make and build the project while the last makes a folder for the solutions. Note that cmake may throw an error and require you to point to the Deal.II library manually.

The code can then be run with the command
```ruby
	./run
```

This will run the code and the solution files will be stored in the "solutions/" folder in vtk format. This can be visualized using a variety of viewers, particularly Paraview and VisIt.

In this problem, we initialize a flat sheet. We then incrementally apply spontaneous curvature and strain to achieve a conical shape. We then hold the spontaneous in-plane strain fixed and incrementally decrease the spontaneous curvature and reverse its sign. This example demonstrates the snap-through instability.