odevo
=====

What is it?
----------

A fast and parallel ordinary differential equation optimizer with compiled systems and easy scripting interfaces. It is particularly suited for stable parameter optimizations since it uses a parallel
version of a stochastic-ranking evolutionary optimizer which should be able to find the best global
parameters consistently.

Why would I need it?
--------------------

`odevo` is particularly helpful if you need to fit parameters of large ordinary differential 
equations systems and other scripting methods are too slow. As such it does not aim to replace 
other great libraries such as (scipy)[www.scipy.org] or (R)[www.r-project.org], but rather gives 
an alternative in case those are to slow. 

How is it faster?
-----------------

Most scripting languages implement optimization as well as the ordinary differential equation
solver via callback functions. Which means, even though the methods are implemented in C, the 
optimized function or ode system are scripted functions which are called from C. Thus, the C
code is calling Pythoon or R functions for instance. This creates significant overhead as those
functions are usually called several thousand times.

`odevo` circumvents the callback overhead by using the only possible solution: the called functions
must be implemented in C(++) as well. However, it is designed in a way that you will have to interact
with C in the lowest amount possible. In fact, for the base application (fitting an ode), the only
thing you have to do is fill in the ode system. This syntax is basically some arithmetic instructions
which are very similar to Python and R. After that odevo will build to executables specific for
your system that allow you to fit, either only the parameters with given initial conditions, or
the parameters and initial conditions to your data set. Calling those programs is than possible
from the command line or minimal helper scripts in the scripting language of your choice.

Requirements
============

You will need the following things installed. **For now there is no support for Windows since
I do not use it and can, therefore, not debug or test it.** 

A C++ and Fortran build environment
---------------------------------

This sounds horrible but it is basically just a compiler for C++ and Fortran. Since it is a vital
part of Linux systems there will be an easy way to install those.

### Debian/Ubuntu

```bash
sudo apt-get install build-essential cmake gfortran
```

### Arch/Manjaro

Arch already comes with a build-system installed, so you only need CMake and Fortran.
 
```bash
sudo pacman -S gcc-fortran cmake
```

MPI
---

In order to use the parallel execution of code we use MPI (Message Passing Interface) which allows
the odevo programs to run on multiple cores or even computer clusters natively. I recommend 
(OpenMPI)[http://www.open-mpi.org/], however, the programs were also tested with 
(MPICH)[http://www.mpich.org/].



Usage
=====

We will first explain the base case (fitting an ode)  
