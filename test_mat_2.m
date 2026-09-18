clear;
clc;


svdTol = 1e-01;
maxK = 50;
A = load_mat_ex_2(); % for test


%--------------------------------------------
%--------------------------------------------
% for printing
markerShapes = {'x', 'p', 'o', 'd'};
colorShapes  = {'b', 'r', 'g', 'm'};
formats = ['d','s', 'h',"q52"];
formats_full = ['fp64','fp32', 'fp16',"q52"];
%--------------------------------------------
%--------------------------------------------


disp('Run teests for Matrix 2')
% Run cases
mat_2_fpDEIM_fp16_randSVD
mat_2_fpDEIM_fp16_svds
mat_2_fpDEIM_fp32_randSVD
mat_2_fpDEIM_fp32_svd
mat_2_fpDEIM_fp64_randSVD
mat_2_fpDEIM_fp64_svd
mat_2_fpDEIM_fp64_svds  
disp('Done')
