function sys_info = pendulum_info(sys_params)

% set default to use stablizing control for path integrals
if(nargin<1)
    use_stable   = true;
    use_unstable = false;
else
    use_stable   = sys_params.use_stable;
    use_unstable = sys_params.use_unstable;
end

%% parameters
m = 1.0;  % Mass (kg)
L = 1.0;  % Length of pendulum (m)
g = 9.81; % Gravity (m/s^2)
c = 0.1;  % Damping coefficient

%% calculate linearization
% x_eqb = [pi 0]


% adding damping = 0.1 to the cart
stable_eqb_pt = 1;
unstable_eqb_pt = -1;

A_damped    = [0 1;+stable_eqb_pt*g/L -c/m*L^2];
A_undamped  = [0 1;-g/L 0];
B           = [0; m*L^2];

[~,D_damped,W_damped]     = eig(A_damped);
[~,D_undamped,W_undamped] = eig(A_undamped);

%% define locally stable system
A = A_damped;

% useling pole placement
if(use_stable)
    k_poles  = place(A,B,[-1.1;-1.2]);
    A_stable = A-B*k_poles;
    sys_info.A_stable = A_stable;
else
    k_poles    = place(A,B,[1.1;1.2]);
    A_unstable = A-B*k_poles;
    sys_info.A_unstable = A_unstable;
end

% using lqr (leads to complex eig vals)
Q     = diag([1 1]);
R     = 1;
k_lqr = lqr(A,B,Q,R);

if(use_stable)
    [~,D,W] = eig(A_stable);
elseif(use_unstable)
    [~,D,W] = eig(A_unstable);
else
    % saddle?
    [~,D,W] = eig(A);
    % realifiy complex eig vals and eig vecs
end
[Wr,Dr] = cdf2rdf(W,D);


%% store system params
sys_info.m  = m;
sys_info.L  = L; 
sys_info.g  = g;
sys_info.c  = c; 

% system setup for baseline control
sys_info.A_undamped = A_undamped;
sys_info.D_undamped = D_undamped;
sys_info.W_undamped = W_undamped;
sys_info.A_damped   = A_damped;
sys_info.D_damped   = D_damped;
sys_info.W_damped   = W_damped;

% linearization params for path integral control
sys_info.A           = A;
sys_info.B           = B;
sys_info.eig_vals    = D;
sys_info.eig_vectors = W;
sys_info.x_eqb       = [pi; 0];

% setup for local control in path integral
sys_info.use_stable     = use_stable;
sys_info.use_unstable   = use_unstable;

% local control 
sys_info.k_poles = k_poles;
sys_info.k_lqr   = k_lqr;

sys_info.id = "pendulum";