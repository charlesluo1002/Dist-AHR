# Distributed Adaptive Huber Regression

## Dist_AHR

Distributed_Huber.R - 
Implements the distributed adaptive huber regression ([Luo, et al., 2022](https://doi.org/10.1016/j.csda.2021.107419)) across local servers for both low and high-dimensional settings. Both asymptotic normal based inference and bootstrap inference are supported.

In low-dimension, it employs Gradient-Descent with Barzilai-Borwein update to optimize the shifted local loss. In high-dimension, iterative local adaptive majorize-minimization algorithm, together with cross validation, is implemented to find the optimal solution.

## Simulation

simulation_low_dim.R - simulates various low-dimensional settings considered in the article to compare the newly proposed dist-AHR against extant methods.

simulation_high_dim.R - simulates various high-dimensional settings considered in the article to compare the newly proposed dist-AHR against extant methods.

## Authors

 Jiyu Luo <jil130@ucsd.edu>, Qiang Sun <qsun@utstat.toronto.edu> and Wen-Xin Zhou <wez243@ucsd.edu>


## Reference
 Luo, J.,  Sun Q. and Zhou, W.-X. (2022). Distributed adaptive Huber regression. *Computational Statistics & Data Analysis*. **169** 107419. [paper](https://doi.org/10.1016/j.csda.2021.107419)
