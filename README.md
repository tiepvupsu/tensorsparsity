# Tensor sparsity 

**Related publication:**
1. Tiep Vu, Lam Nguyen, Calvin Le, Vishal Monga, "[Tensor Sparsity for Classifying Low-Frequency Ultra-Wideband (UWB) SAR Imagery.](http://ieeexplore.ieee.org/document/7944265/)", IEEE Radar Conference 2017. 
2. Tiep Vu, Lam Nguyen, Vishal Monga, "Classifying Multi-channel UWB SAR Imagery via Tensor Sparsity Learning Techniques.", submitted to IEEE Transactions on Aerospace and Electronic Systems.

## Different tensor sparsity constraints 
![alttext](general_SRC.png)

## Main function 
`tensor_sparse_coding.m`
```matlab
function [X] = tensor_sparse_coding(Y, D, Xinit, opts)
% Solving tensor sparse coding problem with different sparsity constraints.
%   using FISTA.
% Y: a 3-D tensor of size (d-n-T), multi-channel signals 
%       (dimension, number of data points, number of channels)
% D: a 3-D tensor of size (d-K-T) dictionary with D 'atoms'
% X: a 3-D tensor of size (K-n-T) starting guess of the solution,
%       if X == [], then it is initialized with all zeros
% opts: a structure, options for the optimization problem
%       opts.eps: a number, tolerance 
%       opts.regul: a string, regularization function:
%           - 'concat': all channels are concatenated to one long long vector (CC)
%           - 'l1': simple l1 (CR)
%           - 'tube': tube sparsity as in the paper (SM)
%           - 'group_tensor': group tensor as in the paper (GT)
%       opts.pos: a boolean, with nonnegativity constraint or not 
% output: 
%       X: a 3-D tensor of size (K-n-T)
```