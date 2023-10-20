# smo-paper
Contains Matlab code to run three numerical examples of a paper submitted to JSMO.

The example files are: ex1.m, ex2.m and ex3.m
These Matlab scripts call CPAS, PENLAB, and FMINSDP to solve an optimization problem.

CPLEX must be installed to run CPAS.
PENLAB can be obtained from: https://web.mat.bham.ac.uk/kocvara/penlab/
and FMINSDP from: https://www.mathworks.com/matlabcentral/fileexchange/43643-fminsdp

See file ex1.m. To run under linux or windows you have to set the path to CPLEX files and execute install.m to prepare PENLAB and FMINSDP.
For example, in my machine:
% % LINUX *************************************************************
%     addpath('/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux')
%     run("/home/marozteg/MATLAB-Drive/fminsdp/install.m")
%     run("/home/marozteg/MATLAB-Drive/PENLABv104/install.m")
% % *******************************************************************
% WINDOWS ***********************************************************
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64')
    run("C:\Users\maroz\MATLAB Drive\fminsdp\install.m")
    run("C:\Users\maroz\MATLAB Drive\PENLABv104\install.m")
% *******************************************************************
