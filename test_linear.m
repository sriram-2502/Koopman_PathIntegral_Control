clc; clear; close all

% defualt plot options
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth', 2) 
set(0,'DefaultAxesFontSize', 20)
set(0,'defaultfigurecolor',[1 1 1])

% add paths
addpath('dynamics')
addpath('baseline_control')
addpath('eigfun_control')
addpath('compute_eigfuns')
addpath('utils')
addpath('animations')


%% setup a random linear system of any dimension
n_states = 3;
x_op = rand(n_states,1);
[~,sys_info] = dynamics_linear(x_op,0);

%% compute path integrals and gradients
phi = setup_path_integrals(x_op, @dynamics_linear);
grad_phi_x_op = compute_gradients(x_op, phi);

%% verify the gradients is the same as left eigenvectors
sys_info.eig_vectors
grad_phi_x_op
error = sys_info.eig_vectors - grad_phi_x_op