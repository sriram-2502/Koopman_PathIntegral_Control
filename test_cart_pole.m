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
sys_params.use_stable   = false; % use forward flow
sys_params.use_unstable = false; % use reverse flow

% add damping to cart
use_damping = true;

%% load dynamics
sys_info = cart_pole_info(sys_params);
dynamics = @dynamics_cart_pole;
n_states = 4; 
n_ctrl   = 1;
x        = sym('x',[n_states;1],'real');
u        = sym('x',[n_ctrl;1],'real');
f        = dynamics_cart_pole(x, u, sys_info);

if(use_damping)
    A = sys_info.A_damped;
    W = sys_info.W_damped;
    D = sys_info.D_damped;
else
    A = sys_info.A_undamped;
    W = sys_info.W_undamped;
    D = sys_info.D_undamped;
end
B = sys_info.B;


%% compute lqr gain
% NOTE: swing up from bottom pos only works when damping is zero (b=0) in
% sys_info setup
Q = diag([200 1000 0 0]);
R  = 0.035;
lqr_params = get_lqr(A,B,Q,R);
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

    % get energy based control
    u = get_swing_up_control(lqr_params, x_op, x_desired);
    use_reverse = false;
    x_next = rk4(@dynamics_cart_pole,dt_sim,x_op',u,use_reverse, sys_info);

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
    x_op = x_next_wrapped;

    % logs 
    Tout  = [Tout;t_sim];
    Xout1 = [Xout1;x_next_wrapped];
    Uout1 = [Uout1;u];

end

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% animate
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

%% state and control plots
figure(22);
% First subplot: x
subplot(2,2,1)
plot(Tout, Xout1(:,1), 'DisplayName', '$x$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('position, $x$', 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on

% Second subplot: x_dot
subplot(2,2,2)
plot(Tout, Xout1(:,2), 'DisplayName', '$\dot{x}$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('velocity, $\dot{x}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid on

% Third subplot: theta
subplot(2,2,3)
plot(Tout, Xout1(:,3), 'DisplayName', '$\theta$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('angle, $\theta$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid

% Fourth subplot: theta_dot
subplot(2,2,4)
plot(Tout, Xout1(:,4), 'DisplayName', '$\dot{\theta}$'); hold on;
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('angular velocity, $\dot{\theta}$', 'Interpreter', 'latex');
box on;
legend('Interpreter', 'latex');
grid