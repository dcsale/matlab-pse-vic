# matlab-VIC #

This is a collection of Matlab code that performs simulation of 3D incompressible flow using a variety of vortex methods.  Currently, a Particle Strength Exchange (PSE) and Vortex-in-Cell (VIC) algorithm are included.  There are some attempts to speedup the code by vectorization, multi-core parallelize, and GPU acceleration.  The PSE method uses direct calculation of Biot-Savart law, so the scaling is O(N^2) and even with these speedups only relatively small problem sizes are possible.  The PSE is a purely particle method, and no resampling (stabilization) of the particles has been finished yet, so the solution can become unstable more quickly than VIC.  The VIC method is based upon high order solution method using FFT solvers and Cartesian background mesh.

Expect many things to be broken and probably accidentally flipped some +'s and -'s at some point ... but have fun anyways!

### What is this repository for? ###

* learning about high order vortex particle methods
* learning about high performance computing
* having fun with Matlab!

### How do I get set up? ###
As of Feb. 17, 2015 ... the code has become a bit unorganized after transitioning between version control methods ... so I hope to improve by providing more clear way to launch specific simulation cases, such as vortex rings, lifting-lines, bluff body flows, etc...

### Contribution guidelines ###

* review the literature included (many PDF papers in documents folder)
* use GIT version control ... commit early and commit often
* other guidelines: ... ?  try to have fun

### Who do I talk to? ###
* contact Danny Clay Sale for questions (dsale@uw.edu)
