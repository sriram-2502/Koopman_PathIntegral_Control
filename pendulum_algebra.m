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
show_animation      = false;
wrap_theta          = true;
show_diagnositcs    = true;

% setup stable or unstable system
sys_params.use_stable   = false; % locallcy stable
sys_params.use_unstable = false; % locally unstable

%% load dynamics
sys_info = pendulum_info(sys_params);
dynamics = @dynamics_pendulum;
n_states = 2; 
n_ctrl   = 1;
x        = sym('x',[n_states;1],'real');
u        = sym('x',[n_ctrl;1],'real');

[dxdt,fx,gx]        = dynamics(x, u, sys_info);
sys_info.dynamics_f = matlabFunction(fx,'vars',x);
sys_info.dynamics_g = matlabFunction(gx,'vars',x);

% get A corresponding to stable or unstable system
if(sys_info.use_stable)
    A = sys_info.A_stable;
else
    A = sys_info.A;
end

B = sys_info.B;
W = sys_info.eig_vectors;
D = sys_info.eig_vals;

% check forward/reverse time path integral
if(all(round(diag(D)) > 0))
    disp('---- using path integrals for unstable system -----')
elseif(all(round(diag(D)) < 0))
    disp('---- using path integrals for stable system -----')
else
    disp('----- using algebra method for saddle system ----')
end

%% compute lqr gain
% baseline for engergy based control
Q_baseline = diag([1 1]);
R_baseline  = 1;
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
x_init      = [pi-0.1 0.0]; 
x_desired   = [pi 0];  
x_eqb       = [pi 0]; 
dt_sim      = 0.01; 
t_start     = 0;
t_end       = 5;
max_iter    = floor(t_end/dt_sim);
x_op1       = x_init;
x_op2       = x_init;
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

% logs
Tout    = t_start;
Xout1   = x_op1; % for baseline controllers and animation
Xout2   = x_op2; % for path integral
Xout1p  = [pi-x_op1(1) x_op1(2)]; % shifted coordinates
Xout2p  = [pi-x_op2(1) x_op2(2)]; % shifted coordinates
Uout1   = []; 
Uout2   = []; 

% diagnostics
integrand_convergence       = []; % plot this to check convergence criteria
eigen_function_estimated    = [];
solution_convergence        = [];

%% start simulation
% show progress bar
w_bar = waitbar(0,'1','Name','running simulation loop...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for t_sim = t_start:dt_sim:t_end
    % udpate progress bar
    waitbar(t_sim/t_end,w_bar,sprintf(string(t_sim)+'/'+string(t_end) +'s'))

    % get current time remaining
    t_span_curr = t_sim:dt_sim:t_end+dt_sim;

    % get energy based control
    u1 = -lqr_params.K_lqr*(x_op1-x_desired)';
    
    % ------ compute eigfn based control ------
    k = 2;
    phi             = compute_path_integrals_algebra(x_op2', dynamics, sys_info, k);
    phi_x_op        = phi.phi_x_op;
    grad_phi_x_op   = compute_gradients(phi);
    P_riccati_curr  = reshape(P_riccati(iter,:),size(A));

    if(show_diagnositcs)
        integrand_convergence = [integrand_convergence;phi.phi_integrand_x_op];
        solution_convergence  = [solution_convergence;phi.phi_sol_conv_x_op];
        eigen_function_estimated = [eigen_function_estimated; phi_x_op];
    end
    
    % u_klqr = compute_control(lqr_params_transformed,P_riccati_curr, phi_x_op, grad_phi_x_op);
    % u_klqr = compute_control_with_riccati(lqr_params,sys_info,phi_x_op,grad_phi_x_op, t_span_curr, x_op2);
    u2     = 0;

    % simulate
    use_reverse = false;
    x_next1 = rk4(dynamics,dt_sim,x_op1',u1,use_reverse,sys_info);
    x_next2 = rk4(dynamics,dt_sim,x_op2',u2,use_reverse,sys_info);

    
    % wrap theta if necessary
    if(wrap_theta)
        theta1 = x_next1(1);
        theta1 = mod(theta1,2*pi);
        x_next1w = [theta1 x_next1(2)];

        theta2 = x_next2(1);
        theta2 = mod(theta2,2*pi);
        x_next2w = [theta2 x_next2(2)];
    else
        x_next1w = x_next1';
        x_next2w = x_next2';
    end

    % shift eqb point for plots
    x_next1p = x_next1w - x_eqb;
    x_next2p = x_next2w - x_eqb;

    % update states
    x_op1 = x_next1w;
    x_op2 = x_next2w;
    iter  = iter + 1;

    % logs 
    Tout   = [Tout;t_sim];
    Xout1  = [Xout1;x_next1w];
    Xout2  = [Xout2;x_next2w];
    Xout1p = [Xout1p;x_next1p];
    Xout2p = [Xout2p;x_next2p];
    Uout1  = [Uout1;u1];
    Uout2  = [Uout2;u2];

end

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% animate
if(show_animation)
    Xanimate = Xout2(:,1);

    hf = figure(11);
    % hf.WindowState = 'maximized';
    skip_rate  = 10;

    % Initialize movieVector array with the correct structure
    numFrames = floor(length(Xanimate) / skip_rate);
    movieVector(1:numFrames) = struct('cdata', [], 'colormap', []);
    
    % Animation loop with frame index
    frameIndex = 1;
    for i = 1:skip_rate:length(Xanimate)
       animate_pendulum(Xanimate(i), sys_info);
       pause(0.1);
    
       % Capture frame and assign it to movieVector
       movieVector(frameIndex) = getframe(hf);
       frameIndex = frameIndex + 1;
    
       hold off
    end
    
    % Remove any empty frames from movieVector (in case loop finished early)
    movieVector = movieVector(1:frameIndex-1);
    
    % Save the movie
    myWriter = VideoWriter('cart_pole', 'MPEG-4');
    myWriter.Quality    = 100;
    myWritter.FrameRate = 180;
    
    % Open the VideoWriter object, write the movie, and class the file
    open(myWriter);
    writeVideo(myWriter, movieVector);
    close(myWriter);
end

%% state and control plots
figure(22);
% First subplot: theta
subplot(2,2,1)
plot(Tout, Xout1p(:,1), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2p(:,1),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('position, $\theta$', 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on

% Second subplot: theta_dot
subplot(2,2,2)
plot(Tout, Xout1p(:,2), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2p(:,2),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('velocity, $\dot \theta$', 'Interpreter', 'latex');
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
    ylabel(['Control ', num2str(i)], 'Interpreter', 'latex');  % Generic label for each state
    legend('Interpreter', 'latex');
    grid on;
end

%%
if(show_diagnositcs)
    figure(44)
    subplot(2,3,1)
    plot(Tout(1:length(integrand_convergence)), integrand_convergence(:,1), 'DisplayName', 'convergence'); hold on;
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('$\exp(-\lambda T) w^\top F_n(s_T(x))$', 'Interpreter', 'latex');

    subplot(2,3,4)
    plot(Tout(1:length(integrand_convergence)), integrand_convergence(:,2), 'DisplayName', 'convergence'); hold on;
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('$\exp(-\lambda T) w^\top F_n(s_T(x))$', 'Interpreter', 'latex');

    subplot(2,3,2)
    eigen_function_linearized = W'*Xout2';
    plot(Tout(1:length(eigen_function_estimated)), eigen_function_estimated(:,1), 'DisplayName', 'estimated'); hold on;
    plot(Tout(1:length(eigen_function_linearized)), eigen_function_linearized(1,:), 'DisplayName', 'linearized'); hold on;
    legend('Interpreter', 'latex');
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('eigen function', 'Interpreter', 'latex');

    subplot(2,3,5)
    eigen_function_estimated_test = eigen_function_estimated;
    eigen_function_estimated_test(abs(eigen_function_estimated_test)>100)=100;
    plot(Tout(1:length(eigen_function_estimated)), eigen_function_estimated(:,2), 'DisplayName', 'estimated'); hold on;
    % plot(Tout(1:length(eigen_function_linearized)), eigen_function_linearized(2,:), 'DisplayName', 'linearized'); hold on;
    legend('Interpreter', 'latex');
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('eigen function', 'Interpreter', 'latex')

    subplot(2,3,3)
    plot(Tout(1:length(solution_convergence)), solution_convergence(:,1), 'DisplayName', 'convergence'); hold on;
    legend('Interpreter', 'latex');
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('solution (X1)', 'Interpreter', 'latex');

    subplot(2,3,6)
    plot(Tout(1:length(solution_convergence)), solution_convergence(:,2), 'DisplayName', 'convergence'); hold on;
    legend('Interpreter', 'latex');
    xlabel('time (s)', 'Interpreter', 'latex');
    ylabel('solution (X2)', 'Interpreter', 'latex')

end
