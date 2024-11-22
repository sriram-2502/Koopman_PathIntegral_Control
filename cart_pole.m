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

% setup stable or unstable system
sys_params.use_stable   = true; % use forward flow
sys_params.use_unstable = false; % use reverse flow

% fix seed
rng(15)

%% setup a random linear system of any dimension
dynamics        = @dynamics_cart_pole;
n_states        = 4; 
n_ctrl          = 1; % TODO: check wth n_ctrl > 1
x_op            = rand(n_states,1);
x               = sym('x',[n_states;1],'real');
u               = sym('x',[n_ctrl;1],'real');
[f,sys_info]    = dynamics(x, u, sys_params);
B               = sys_info.B;
W               = sys_info.eig_vectors;
D               = sys_info.eig_vals; %TODO: check if order matters

% check forward/reverse time path integral
if(all(ceil(diag(D)) <= 0))
    disp('---- using forward time path integrals -----')
elseif(all(ceil(diag(D)) > 0))
    disp('---- using reverse time path integrals -----')
end

% get A corresponding to stable or unstable system
if(sys_info.use_stable)
    A = sys_info.A_stable;
elseif(sys_info.use_unstable)
    A = sys_info.A_unstable;
end


%% compute lqr gain
Q = diag([200 1000 0 0]);
R  = 0.035;
lqr_params = get_lqr(A,B,Q,R);

% solving lqr in koopman coordinates
A_transformed           = D;
B_transformed           = W'*B;
Q_transformed           = inv(W)*Q*inv(W');
lqr_params_transformed  = get_lqr(A_transformed,B_transformed,Q_transformed,R);

%% simulation loop
x_init      = [0.0 pi-0.1 0.0 0.0]; 
x_desired   = [0.0 pi 0.0 0.0];  
x_eqb       = [0.0 pi 0.0 0.0]; 
dt_sim      = 0.01; 
t_start     = 0;
t_end       = 1;
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
Xout1   = [x_op(1) pi-x_op(2) x_op(3) x_op(4)]; % shifted coordinates for path integrals
Xout2   = x_op; % for baseline controllers and animation
Uout1   = []; 

% show progress bar
w_bar = waitbar(0,'1','Name','running simulation loop...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for t_sim = t_start:dt_sim:t_end

    % udpate progress bar
    iter = ceil(t_sim/t_end)+1;
    waitbar(t_sim/t_end,w_bar,sprintf(string(t_sim)+'/'+string(t_end) +'s'))
    
    % ------ compute eigfn based control ------
    phi             = compute_path_integrals(x_op1', dynamics, sys_params);
    phi_x_op        = phi.phi_x_op;
    grad_phi_x_op   = compute_gradients(phi);
    P_riccati_curr  = reshape(P_riccati(iter,:),size(A));
    u1              = compute_control(lqr_params_transformed,P_riccati_curr, phi_x_op, grad_phi_x_op);

    % get energy based control
    u2 = get_swing_up_control(@dynamics_cart_pole, lqr_params, x_op, x_desired);

    % simulate the system
    use_reverse = false; % do forward simulation for control loop
    x_next1     = euler(dynamics,dt_sim,x_op1',u1,use_reverse,sys_params);
    x_next2     = euler(dynamics,dt_sim,x_op2',u2,use_reverse,sys_params);

    % wrap theta if necessary
    if(wrap_theta)
        theta = x_next2(2);
        if(abs(theta)<1e-3)
            theta = 0;
        end
        theta = mod(theta,2*pi);
        x_next_wrapped = [x_next2(1) theta x_next2(3) x_next2(4)];
    else
        x_next_wrapped = x_next2';
    end

    % shift eqb point to unstable point
    x_next_saddle = x_next1' - x_eqb;
    
    % ------ update states ------
    x_op1 = x_next1';
    x_op2 = x_next2';
    iter  = iter + 1;

    % logs 
    Tout  = [Tout;t_sim];
    Xout1 = [Xout1;x_next_saddle];
    Xout2 = [Xout2;x_next_wrapped'];
    Uout1 = [Uout1;u1];
    Uout2 = [Uout2;u2];

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