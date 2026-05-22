The folder contains an implementation of the algorithm and the examples submitted to EMSOFT 2016. Each example comes with a demo file start with main_*.

To run the source code, two libraries needs to be installed on MATLAB: yalmip and SDPT3 (or Sedumi), which can be found using the following link:
http://users.isy.liu.se/johanl/yalmip/
http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
http://coral.ise.lehigh.edu/software/

If you see errors like "Undefined function 'sdpsettings' for input arguments of type 'char'.", please go back and check the required libraries.
Note different optimization toolbox make cause slightly different results.

The current implementation now contains an MATLAB ode solver for simulations.
For guaranteed simulate, one needs to install library CAPD to generate traces, which will replace the traces computed by the built in MATLAB solver. The CAPD library can be found at: http://capd.ii.uj.edu.pl/
