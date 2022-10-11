# LeapingLCEs
Modeling deformations of a LCE sheet

## Installation
This code relies on Deal.II to do finite element calculations. For detailed instructions and documentation on this software, refer to https://www.dealii.org/

After Deal.II is installed, navigate to the source directory (on the level of the license and README files) and enter the following commands:

mkdir build && cd build/
cmake ..
make
mkdir solutions

The first three commands make and build the project while the last makes a folder for the solutions. Note that cmake may throw an error and require you to point to the Deal.II library manually.