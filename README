Installation
------------

Required libraries: gsl, gslblas, pthread

On Linux/Unix run

 ./configure
 make; make install

On Mac OS, the location of the required gsl, gslblas and pthread
libraries may need to be specified:

 ./configure LDFLAGS="-L/opt/local/lib" CPPFLAGS="-I/opt/local/include"
 make; make install

The binary 'gaprec' will be installed in /usr/local/bin unless a
different prefix is provided to configure. (See INSTALL.)

HGAPREC: Hierarchical Gamma Poisson factorization based recommendation tool
----------------------------------------------------------------------------

**hgaprec** [OPTIONS]

   -dir <string>    path to dataset directory with 3 files:
   		    train.tsv, test.tsv, validation.tsv
		    (for examples, see example/movielens-1m)
 
   -m <int>	  number of items
   -n <int>	  number of users
   -k <int>	  number of factors
   
   -rfreq <int>	  assess convergence and compute other stats 
   		  <int> number of iterations
		  default: 10

   -a
   -b		  set hyperparameters
   -c		  default: a = b = c = d = 0.3
   -d
   
   
   -hier	  learn the hierarchical model with Gamma priors
                  on user and item scale parameters

   -bias	  use user and item bias terms
   
   -binary-data	  treat observed data as binary
   		  (if rating > 0 then rating is treated as 1)
   		  
   -gen-ranking	  generate ranking file to use in precision 
   		  computation; see example		  

   -msr		  write out ranking file assuming the test file
                  is based on leave-one-out, i.e., leaving one
                  item out for each user

Example
--------

1. Input data

