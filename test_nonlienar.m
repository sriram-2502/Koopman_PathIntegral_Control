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
show_diagnositcs    = true;

% use linearization for solving riccati in transformed coordinates
sys_params.use_linear_riccati = true; 

%% load dynamics
sys_info = nonlinear_sys_info(sys_params);
dynamics = @dynamics_nonlinear;
n_states = sys_info.state_dim; 
n_ctrl   = sys_info.ctrl_dim;
x        = sym('x',[n_states;1],'real');
u        = sym('x',[n_ctrl;1],'real');
dx_dxt   = dynamics(x, u, sys_info);

A = sys_info.A;
B = sys_info.B;
W = sys_info.eig_vectors;
D = sys_info.eig_vals;

% check forward/reverse time path integral
if(all(round(diag(D)) >= 0))
    disp('---- using path integrals for unstable system -----')
elseif(all(round(diag(D)) < 0))
    disp('---- using path integrals for stable system -----')
else
    disp('---- using algebra for saddle system ----')
end

%% compute lqr gain
Q_baseline = eye(size(n_states));
R_baseline  = ones(size(n_ctrl));
lqr_params_baseline = get_lqr(sys_info.A,B,Q_baseline,R_baseline);

% for klqr - path integral based control
Q           = Q_baseline;
R           = R_baseline;
lqr_params  = get_lqr(A,B,Q,R);

% solving lqr in koopman coordinates
A_transformed           = D;
B_transformed           = W'*B;
Q_transformed           = inv(W)*Q*inv(W');
lqr_params_transformed  = get_lqr(A_transformed,B_transformed,Q_transformed,R);

%% simulation loop
x_init      = 10*rand(n_states,1);  
dt_sim      = 0.01; 
t_start     = 0;
t_end       = 5;
max_iter    = floor(t_end/dt_sim);
x_op1       = x_init';
x_op2       = x_init';
t_span      = t_start:dt_sim:t_end;
iter        = 1;

% solve differential riccati
[t_riccati,P_riccati] = compute_riccati(lqr_params_transformed, t_span);
P_riccati_curr = reshape(P_riccati(1,:),size(A));

% check P matrices are the same
if(show_diagnositcs)
    % check P in transformed coordinates for finite vs inf horizon
    disp('check P in transformed coordinates for finite vs inf horizon')
    P_infinite = lqr_params_transformed.P_lqr;
    P_finite   = P_riccati_curr;
    riccati_error = P_infinite-P_finite
    
    % check P in both transformed and orig coordintes for inf horizon
    disp('Transform P to orig coordinates and check error for inf horizon')
    [K,P,e] = lqr(A,B,Q,R);
    P_error = P - W*P_infinite*W'
end

%% start simulation
% logs
Tout    = t_start;
Xout1   = x_op1; % for baseline controllers and animation
Xout2   = x_op2; % 
Uout1   = []; 
Uout2   = []; 
convergence = []; % plot this to check convergence criteria

% show progress bar
w_bar = waitbar(0,'1','Name','running simulation loop...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for t_sim = t_start:dt_sim:t_end
    % udpate progress bar
    waitbar(t_sim/t_end,w_bar,sprintf(string(t_sim)+'/'+string(t_end) +'s'))

    % get current time remaining
    t_span_curr = t_sim:dt_sim:t_end+dt_sim;

    % get baseline lqr control
    u1 = -lqr_params_baseline.K_lqr*x_op1';
    
    % ------ compute eigfn based control ------
    phi             = compute_path_integrals(x_op2', dynamics, sys_info);
    phi_x_op        = phi.phi_x_op;
    grad_phi_x_op   = compute_gradients(phi);
    P_riccati_curr  = reshape(P_riccati(iter,:),size(A));
    u2              = compute_control(lqr_params_transformed,P_riccati_curr, phi_x_op, grad_phi_x_op);
    
    if(show_diagnositcs)
        convergence = [convergence;phi.phi_integrand_x_op];
    end

    % simulate
    use_reverse = false;
    x_next1 = rk4(dynamics,dt_sim,x_op1',u1,use_reverse,sys_info);
    x_next2 = rk4(dynamics,dt_sim,x_op2',u2,use_reverse,sys_info);


    % update states
    x_op1 = x_next1';
    x_op2 = x_next2';
    iter  = iter + 1;

    % logs 
    Tout   = [Tout;t_sim];
    Xout1  = [Xout1;x_next1'];
    Xout2  = [Xout2;x_next2'];
    Uout1  = [Uout1;u1];
    Uout2  = [Uout2;u2];

end

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);


%% state and control plots
figure(22);
% First subplot: x
subplot(2,2,1)
plot(Tout, Xout1(:,1), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2(:,1),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('position, $x$', 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on

% Second subplot: theta
subplot(2,2,2)
plot(Tout, Xout1(:,2), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2(:,2),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('angle, $\theta$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid on

figure(33);
% Loop over each control and create subplots
for i = 1:n_ctrl
    % Create a subplot grid
    subplot(n_ctrl, 1, i);
    
    % Plot the state data for Xout1 and Xout2
    plot(Tout(1:length(Uout1)), Uout1(:, i),  '-*', 'DisplayName', 'baseline'); hold on;
    plot(Tout(1:length(Uout2)), Uout2(:, i), '-*',  'DisplayName', 'klqr'); hold on;
    
    % Set labels and legend
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel(['Control ', num2str(i)], 'Interpreter', 'latex'); 
    legend('Interpreter', 'latex');
    grid on;
end

if(show_diagnositcs)
    figure(44)
    subplot(2,1,1)
    plot(Tout(1:length(convergence)), convergence(:,1), 'DisplayName', 'baseline'); hold on;
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('$\exp(-\lambda T) w^\top F_n(s_T(x))$', 'Interpreter', 'latex');

    subplot(2,1,2)
    plot(Tout(1:length(convergence)), convergence(:,2), 'DisplayName', 'baseline'); hold on;
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('$\exp(-\lambda T) w^\top F_n(s_T(x))$', 'Interpreter', 'latex');
end
