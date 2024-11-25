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
addpath('solve_riccati')
addpath('utils')
addpath('animations')

% setup params
show_animation      = true;
wrap_theta          = true;
show_diagnositcs    = true;

% setup stable or unstable system
sys_params.use_stable   = true; % use forward flow
sys_params.use_unstable = false; % use reverse flow

%% load dynamics
dynamics        = @dynamics_cart_pole;
n_states        = 4; 
n_ctrl          = 1;
x               = sym('x',[n_states;1],'real');
u               = sym('x',[n_ctrl;1],'real');
[f,sys_info]    = dynamics_cart_pole(x, u, sys_params);
B               = sys_info.B;
W               = sys_info.eig_vectors;
D               = sys_info.eig_vals;

% check forward/reverse time path integral
if(all(ceil(diag(D)) <= 0))
    disp('---- using forward time path integrals -----')
elseif(all(ceil(diag(D)) > 0))
    disp('---- using reverse time path integrals -----')
end

% get A corresponding to stable or unstable system
if(sys_info.use_stable)
    A = sys_info.A_stable;
else
    A = sys_info.A;
end

%% compute lqr gain
% baseline for engergy based control
Q_baseline = diag([200 1000 0 0]);
R_baseline  = 0.035;
lqr_params_baseline  = get_lqr(sys_info.A,B,Q_baseline,R_baseline);

% for klqr - path integral based control
Q           = diag([200 1000 0 0]);
R           = 0.035;
lqr_params  = get_lqr(A,B,Q,R);

% solving lqr in koopman coordinates
A_transformed           = D;
B_transformed           = W'*B;
Q_transformed           = inv(W)*Q*inv(W');
lqr_params_transformed  = get_lqr(A_transformed,B_transformed,Q_transformed,R);

%% simulation loop
x_init      = [0.0 0.1 0.0 0.0]; 
x_desired   = [0.0 pi 0.0 0.0];  
x_eqb       = [0.0 pi 0.0 0.0]; 
dt_sim      = 0.01; 
t_start     = 0;
t_end       = 10;
max_iter    = floor(t_end/dt_sim);
x_op1       = x_init;
x_op2       = [x_init(1) pi-x_init(2) x_init(3) x_init(4)];
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
Xout2   = x_op2; % 
Xout1p  = [x_op1(1) pi-x_op1(2) x_op1(3) x_op1(4)]; % shifted coordinates
Xout2p  = [x_op1(1) pi-x_op1(2) x_op1(3) x_op1(4)]; % shifted coordinates
Uout1   = []; 
Uout2   = []; 

% show progress bar
w_bar = waitbar(0,'1','Name','running simulation loop...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for t_sim = t_start:dt_sim:t_end

    % udpate progress bar
    iter = ceil(t_sim/t_end)+1;
    waitbar(t_sim/t_end,w_bar,sprintf(string(t_sim)+'/'+string(t_end) +'s'))

    % get energy based control
    u1 = get_swing_up_control(@dynamics_cart_pole, lqr_params_baseline, x_op1, x_desired);
    
    % ------ compute eigfn based control ------
    phi             = compute_path_integrals(x_op2', dynamics, sys_params);
    phi_x_op        = phi.phi_x_op;
    grad_phi_x_op   = compute_gradients(phi);
    P_riccati_curr  = reshape(P_riccati(iter,:),size(A));
    u_volt          = compute_control(lqr_params_transformed,P_riccati_curr, phi_x_op, grad_phi_x_op);
    u_volt          = -sys_info.k_poles*x_op2' + u_volt;
    u_volt          = saturate_fun(u_volt,12,-12);
    x_dot           = x_op2(3);
    u_klqr          = volts_to_force(x_dot,u_volt);
    u2              = u_klqr;

    % simulate
    use_reverse = false;
    x_next1 = rk4(dynamics,dt_sim,x_op1',u1,use_reverse,sys_params);
    x_next2 = rk4(dynamics,dt_sim,x_op2',u2,use_reverse,sys_params);

    % wrap theta if necessary
    if(wrap_theta)
        theta1 = x_next1(2);
        if(abs(theta1)<1e-3)
            theta1 = 0;
        end
        theta1 = mod(theta1,2*pi);
        x_next1w = [x_next1(1) theta1 x_next1(3) x_next1(4)];

        theta2 = x_next2(2);
        if(abs(theta2)<1e-3)
            theta2 = 0;
        end
        theta2 = mod(theta2,2*pi);
        x_next2w = [x_next2(1) theta2 x_next2(3) x_next2(4)];
    else
        x_next1w = x_next1';
        x_next2w = x_next2';
    end

    % shift eqb point to unstable point
    x_next1p = x_next1' - x_eqb;
    x_next2p = x_next2' - x_eqb;
    
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
    Xanimate = Xout2;
    % Xanimate(:,1) = 0;
    Xanimate(:,2) = Xanimate(:,2) + pi;
    hf = figure(11);
    % hf.WindowState = 'maximized';
    skip_rate  = 10;
    
    % Initialize movieVector array with the correct structure
    numFrames = floor(length(Xanimate) / skip_rate);
    movieVector(1:numFrames) = struct('cdata', [], 'colormap', []);
    
    % Animation loop with frame index
    frameIndex = 1;
    for i = 1:skip_rate:length(Xanimate)
       animate_cart_pole(Xanimate(i,1),Xanimate(i,2));
       pause(0.01);
    
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
% First subplot: x
subplot(2,2,1)
plot(Tout, Xout1p(:,1), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2p(:,1),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('position, $x$', 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on

% Second subplot: theta
subplot(2,2,2)
plot(Tout, Xout1p(:,2), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2p(:,2),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('angle, $\theta$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid on

% Third subplot: x_dot
subplot(2,2,3)
plot(Tout, Xout1p(:,3), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2p(:,3),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('velocity, $\dot{x}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid

% Fourth subplot: theta_dot
subplot(2,2,4)
plot(Tout, Xout1p(:,4), 'DisplayName', 'baseline'); hold on;
plot(Tout, Xout2p(:,4),'-', 'DisplayName', 'klqr'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('angular velocity, $\dot{\theta}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid


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