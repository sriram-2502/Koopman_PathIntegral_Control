clc; clear; close all
addpath('dynamics')
addpath('baseline_control')
addpath('compute_eigfun')
addpath('solve_riccati')
addpath('utils')

%% load dynamics
n_states = 4; 
n_ctrl   = 1;
x        = sym('x',[n_states;1],'real');
u        = sym('x',[n_ctrl;1],'real');
sys_info = dynamics_cart_pole(x, u);

% utils function
sat = @(x, x_max, x_min) min( x_max, max(x_min,x) ); % Saturation Function

%% compute lqr gain
Q = diag([200 1000 0 0]);
R  = 0.035;
lqr_params = get_lqr(sys_info.A,sys_info.B,Q,R);
K_lqr = lqr_params.K_lqr;

%% simulation loop
x_init      = [0.0 0.1 0.0 0.0]; 
x_desired   = [0.0 pi 0.0 0.0];            
dt_sim      = 0.01; 
t_start     = 0;
t_end       = 10;
max_iter    = floor(t_end/dt_sim);
x_op        = x_init;

% logs
Tout    = t_start;
Xout1   = x_op; 
Uout1   = [];

% show progress bar
w_bar = waitbar(0,'1','Name','running simulation loop...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for t_sim = t_start:dt_sim:t_end

    % udpate progress bar
    waitbar(t_sim/t_end,w_bar,sprintf(string(t_sim)+'/'+string(t_end) +'s'))

    u = 0;

    x_next = rk4(@dynamics_cart_pole,dt_sim,x_op,u);
    theta = x_next(2);

    % take theta as (360- theta), when theta < 0, otherwise theta = theta
    % only for LQR control 
    if theta < 0
       theta = 2*pi-abs(theta);
       x_next_wrapped = [x_next(1) theta x_next(3) x_next(4)];
    else
       x_next_wrapped = x_next';
    end
    
    % update states
    x_op = x_next';

    % logs 
    Tout  = [Tout;t_sim];
    Xout1 = [Xout1;x_next_wrapped];
    Uout1 = [Uout1;u];

end

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% state and control plots
figure(11);
% First subplot: x
subplot(2,2,1)
plot(Tout, Xout1(:,1), 'DisplayName', '$x$'); hold on;
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State: $x$', 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on

% Second subplot: x_dot
subplot(2,2,2)
plot(Tout, Xout1(:,2), 'DisplayName', '$\dot{x}$'); hold on;
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State: $\dot{x}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid on

% Third subplot: theta
subplot(2,2,3)
plot(Tout, Xout1(:,3), 'DisplayName', '$\theta$'); hold on;
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State: $\theta$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid

% Fourth subplot: theta_dot
subplot(2,2,4)
plot(Tout, Xout1(:,4), 'DisplayName', '$\dot{\theta}$'); hold on;
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State: $\dot{\theta}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid