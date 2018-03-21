# bnbmt
multi-threaded bnb method 

Building 
========

First edit all.inc, then run

make dep
make

Running
=======

./<name>.exe [np tres maxsteps]

optional parameters:

np -- number of virtual processors
tres -- the minimal number of iterations to start parallelization (below this value - use one thread)
maxsteps -- maximal total number of iterations