# SRHC Implementation
We have implemented the algorithm presented in the paper as a part of our course course project


![](Figures/State_norm_norm0_40.png)
## Acknowledgement :
This work was developed and proven by the authors of the paper : \
[Stochastic receding horizon control with output feedback and bounded controls](https://www.sciencedirect.com/science/article/abs/pii/S0005109811004882)
### Authors : 
* Peter Hokayem
* Eugenio Cinquemani
* Debasish Chatterjee
* Federico Ramponi
* John Lygeros

## Instructions for running the code : 
* Make sure YALMIP and sdpt3 packages are installed in MATLAB to run the code 
* Phi_func applies piecewise linear saturation to bound values between $-\varphi_{max}$ and $\varphi_{max}$
* Monte carlo simulation file simulates $10^5$ trajectories to compute $\Lambda_t$ values 
* Reachability_matrix computes $\mathcal{R}_k(A,B) = [ A^{k-1}B\ \dots\ AB\ B ]$
* main_lqg performs clipped lqg on system described in it
* main_new and main_new_2 performs algorithm 1 on examples 1 and 2 as described in the paper
* main_new_alg2 performs algorithm 2 on example 1 
