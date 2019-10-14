% This comparison use the Tensor Toolbox supplemental software
% To get this toolbox, visit :http://www.tensortoolbox.org/. 

% Change Rotation Algorithm to do the sort of norm at the same time in paper:
% <On Parallel Implementation of the One-sided Jacobi Algorithm>

% Ring Jacobi array and its proof in paper:
% <A Parallel Ring Ordering Algorithm for EfficientOne-Sided Jacobi SVD
% Computations>

% Iteration way and its proof in paper:
% <Iterative Superlinear-Convergence SVD Beamforming Algorithm and VLSI
% Architecture for MIMO-OFDM Systems>

% Another way I'm still trying to dieive leading orthogonal vector without 
% implementing full SVD in paper:
% <Improving the speed of multi-wau algorithms in Tucker3>

% Reference Paper:
% <Tucker Decomposition for FPGA>

%% initialization
info=create_problem('Type','Tucker','Size',[128,128,128],'Num_Factors',[32,32,32]);
N=3;
order=32;
init=cell(N,1);
for i=1:3
    temp=eye(128);
    init{i}=temp(:,1:order);
end

%% HOOI test without noise


% reference
T=tucker_als(info.Soln,[32,32,32]);

% my work
T1=HOOI(info.Soln,order,init,1);

%% HOOI test with noised data

% reference
NT=tucker_als(info.Data,[32,32,32]);
ALLT=full(NT);

% my work
NT1=HOOI(info.Data,order,init,1);
ALLT1=full(NT1);

%% HOOI using ite_svd with noised data

NT2=HOOI(info.Data,order,init,0);
ALLT2=full(NT2);