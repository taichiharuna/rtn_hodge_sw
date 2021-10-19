# rtn_hodge_sw
Repository for Haruna, T., Fujiki, Y., 2016. Hodge decomposition of information flow on small-world networks. Frontiers in Neural Circuits 10, 77.
The final authenticated version of the paper is available online at:https://www.frontiersin.org/articles/10.3389/fncir.2016.00077/full

The followings are the source code and data used in Figure 2 of the paper.

We use GCC (https://gcc.gnu.org/) and GSL (https://www.gnu.org/software/gsl/). 

## rtn_hodge_sw_demo.c
The code for simulating random threshold networks (RTNs) on the WS model and calculating the Hodege decomposition of the information flow associated with the RTN dynamics. This code can be used to reproduce the other results in the paper by chaging the parameters. 

## text files
Outputs of the rtn_hodge_sw_demo.c. In particular, 

### edge flow.txt
contains the information flow in the matrix form. Its components in the Hodge decomposition are "gradientflow.txt", "curlflow.txt" and "harmonicflow.txt". 

### rtn_hodge_sw_trial_n8_k2.txt
contains the l^2-norms of the edge, gradient, harmonic and curl flows in the columns 4, 5, 6 and 7, respectively.
