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
show_animation = true;
wrap_theta = true;
show_diagnositcs = true;

%% load dynamics
n_states        = 4; 
n_ctrl          = 1;
x               = sym('x',[n_states;1],'real');
u1              = sym('x',[n_ctrl;1],'real');
[f,sys_info]    = dynamics_cart_pole(x, u1);
A               = sys_info.A;
B               = sys_info.B;
W               = sys_info.eig_vectors;
D               = sys_info.eig_vals;

%% compute lqr gain
Q = diag([200 1000 0 0]);
R  = 0.035;
lqr_params = get_lqr(sys_info.A,sys_info.B,Q,R);

%% simulation loop
x_init      = [0.0 pi-0.1 0.0 0.0]; 
x_desired   = [0.0 pi 0.0 0.0];  
x_eqb       = [0.0 pi 0.0 0.0]; 
dt_sim      = 0.01; 
t_start     = 0;
t_end       = 1;
max_iter    = floor(t_end/dt_sim);
x_op        = x_init;
t_span      = t_start:dt_sim:t_end;

% solve differential riccati
Q = eye(size(A));
R = 0.1;
lqr_params_transformed = get_lqr_transformed(D,W'*B,Q,R);
[t_riccati,P_riccati] = compute_riccati(sys_info, lqr_params_transformed, t_span);
P_riccate_cur = reshape(P_riccati(1,:),size(A));

% check P matrices are the same
if(show_diagnositcs)
    % solving lqr in koopman coordinates
    [K_lqr,P_lqr,e] = lqr(D,W'*B,lqr_params_transformed.Q,lqr_params_transformed.R);
    disp('check if P_lqr matches with P_riccati')
    P_lqr
    P_riccate_cur
end

% logs
Tout    = t_start;
Xout1   = x_op; % for baseline controllers and animation
Xout2   = [x_op(1) pi-x_op(2) x_op(3) x_op(4)]; % shifted coordinates
Uout1   = []; 

% show progress bar
w_bar = waitbar(0,'1','Name','running simulation loop...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for t_sim = t_start:dt_sim:t_end

    % udpate progress bar
    iter = ceil(t_sim/t_end)+1;
    waitbar(t_sim/t_end,w_bar,sprintf(string(t_sim)+'/'+string(t_end) +'s'))

    % get energy based control
    u1 = get_swing_up_control(@dynamics_cart_pole, lqr_params, x_op, x_desired);
    
    % compute eigfn based control
    phi = setup_path_integrals(x_op, @dynamics_cart_pole);
    grad_phi = compute_gradients(x_op, phi);
    P_riccati_curr = reshape(P_riccati(iter,:),size(A));
    u2 = compute_control(sys_info, lqr_params_transformed, P_riccati_curr, phi, grad_phi);

    % simulate using rk4
    x_next = rk4(@dynamics_cart_pole,dt_sim,x_op,0);

    % wrap theta if necessary
    if(wrap_theta)
        theta = x_next(2);
        if(abs(theta)<1e-3)
            theta = 0;
        end
        theta = mod(theta,2*pi);
        x_next_wrapped = [x_next(1) theta x_next(3) x_next(4)];
    else
        x_next_wrapped = x_next';
    end

    % shift eqb point to unstable point
    x_next_saddle = x_next' - x_eqb;
    
    % update states
    x_op = x_next';

    % logs 
    Tout  = [Tout;t_sim];
    Xout1 = [Xout1;x_next_wrapped];
    Xout2 = [Xout2;x_next_saddle];
    Uout1 = [Uout1;u1];

end

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% animate
if(show_animation)
    hf = figure(11);
    % hf.WindowState = 'maximized';
    skip_rate  = 10;
    
    % Initialize movieVector array with the correct structure
    numFrames = floor(length(Xout1) / skip_rate);
    movieVector(1:numFrames) = struct('cdata', [], 'colormap', []);
    
    % Animation loop with frame index
    frameIndex = 1;
    for i = 1:skip_rate:length(Xout1)
       animate_cart_pole(Xout1(i,1),Xout1(i,2));
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
Xout = Xout2;
figure(22);
% First subplot: x
subplot(2,2,1)
plot(Tout, Xout(:,1), 'DisplayName', '$x$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('position, $x$', 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on

% Second subplot: theta
subplot(2,2,2)
plot(Tout, Xout(:,2), 'DisplayName', '$\theta$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('angle, $\theta$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid on

% Third subplot: x_dot
subplot(2,2,3)
plot(Tout, Xout(:,3), 'DisplayName', '$\dot{x}$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('velocity, $\dot{x}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid

% Fourth subplot: theta_dot
subplot(2,2,4)
plot(Tout, Xout(:,4), 'DisplayName', '$\dot{\theta}$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('angular velocity, $\dot{\theta}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid