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

% setup params
show_animation      = true;
wrap_theta          = true;
show_diagnositcs    = true;

% fix seed
rng(15)

%% setup a random linear system of any dimension
dynamics        = @dynamics_linear;
n_states        = 4; 
n_ctrl          = 1;
x_op1           = rand(n_states,1);
x               = sym('x',[n_states;1],'real');
u               = sym('x',[n_ctrl;1],'real');
[f,sys_info]    = dynamics(x, u);
A               = sys_info.A_stable;
B               = sys_info.B;
W               = sys_info.eig_vectors;
D               = sys_info.eig_vals;

if(all(diag(D) < 0))
    disp('---- using forward time path integrals -----')
elseif(all(diag(D) > 0))
    disp('---- using reverse time path integrals -----')
end

%% compute path integrals and gradients
phi             = setup_path_integrals(x_op1, dynamics);
grad_phi_x_op   = compute_gradients(x_op1, phi);

%% verify the gradients is the same as left eigenvectors
disp('check eigenfunction gradient matches with left eigenvector of A')
gradient_error = sys_info.eig_vectors - grad_phi_x_op

%% compute lqr gain
Q             = eye(n_states,n_states);
Q_transformed = inv(W')*Q*W';
R             = ones(n_ctrl);

% solving lqr in koopman coordinates
A_transformed = D;
B_transformed = W'*B;
lqr_params    = get_lqr(A_transformed,B_transformed,Q_transformed,R);

%% simulation loop
x_init      = 10*rand(n_states,1);
dt_sim      = 0.1; 
t_start     = 0;
t_end       = 5;
max_iter    = floor(t_end/dt_sim);
x_op1       = x_init';
x_op2       = x_init';
t_span      = t_start:dt_sim:t_end;

% solve differential riccati
[t_riccati,P_riccati] = compute_riccati(sys_info, lqr_params, t_span);
lqr_params.P_riccati_curr = reshape(P_riccati(1,:),size(A));

% check P matrices are the same
if(show_diagnositcs)
    disp('check if P_lqr matches with P_riccati')
    riccati_error = lqr_params.P_lqr-lqr_params.P_riccati_curr
end

% logs
Tout    = t_start;
Xout1   = x_op1;
Xout2   = x_op2;
Uout1   = []; 
Uout2   = []; 

if(show_diagnositcs)
    phi_error = [];
end

% show progress bar
p_bar = waitbar(0,'1','Name','running simulation loop...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for t_sim = t_start:dt_sim:t_end

    % udpate progress bar
    iter = ceil(t_sim/t_end)+1;
    waitbar(t_sim/t_end,p_bar,sprintf(string(t_sim)+'/'+string(t_end) +'s'))
    
    % compute eigfn based control
    phi                         = setup_path_integrals(x_op1', dynamics);
    phi_x_op                    = phi.phi_x_op;
    grad_phi_x_op               = compute_gradients(x_op1', phi);
    lqr_params.P_riccati_curr   = reshape(P_riccati(iter,:),size(A));
    u1                          = compute_control(sys_info, lqr_params, phi_x_op, grad_phi_x_op);

    % get baseline lqr control
    K  = lqr(A,B,Q,R);
    u2 = K*x_op2';

    % simulate the system
    x_next1 = rk4(dynamics,dt_sim,x_op1',u1);
    x_next2 = rk4(dynamics,dt_sim,x_op2',u2);
  
    % update states
    x_op1 = x_next1';
    x_op2 = x_next2';

    % logs 
    Tout  = [Tout;t_sim];
    Xout1 = [Xout1;x_op1];
    Xout2 = [Xout2;x_op2];
    Uout1 = [Uout1;u1];
    Uout2 = [Uout2;u2];

    if(show_diagnositcs)
        phi_error = phi_x_op - (W'*x_op2')';
        phi_error = [phi_error;phi_error];
    end

end

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% state and control plots
figure(22);
% Loop over each state and create subplots
for i = 1:n_states
    % Create a subplot grid
    subplot(n_states, 1, i);
    
    % Plot the state data for Xout1 and Xout2
    plot(Tout, Xout1(:, i), '-*',  'DisplayName', 'klqr'); hold on;
    plot(Tout, Xout2(:, i), '--*', 'DisplayName', 'lqr'); hold on;
    
    % Set labels and legend
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel(['State ', num2str(i)], 'Interpreter', 'latex');  % Generic label for each state
    legend('Interpreter', 'latex');
    grid on;
end

figure(33);
% Loop over each control and create subplots
for i = 1:n_ctrl
    % Create a subplot grid
    subplot(n_ctrl, 1, i);
    
    % Plot the state data for Xout1 and Xout2
    plot(Tout(1:length(Uout1)), Uout1(:, i),  '-*', 'DisplayName', 'klqr'); hold on;
    plot(Tout(1:length(Uout2)), Uout2(:, i), '--*',  'DisplayName', 'lqr'); hold on;
    
    % Set labels and legend
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel(['Control ', num2str(i)], 'Interpreter', 'latex');  % Generic label for each state
    legend('Interpreter', 'latex');
    grid on;
end