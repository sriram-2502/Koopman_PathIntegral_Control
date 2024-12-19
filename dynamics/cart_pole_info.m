function sys_info = cart_pole_info(sys_params)

% set default to use stablizing control for path integrals
if(nargin<1)
    use_stable   = true;
    use_unstable = false;
else
    use_stable   = sys_params.use_stable;
    use_unstable = sys_params.use_unstable;
end

%% parameters
r   = 0.006;       % radius of motor shaft
M   = 0.135;       % Mass of cart
m   = 0.1;         % mass of pendulum
I   = 0.0007176;   % MOI of Pendulum
l   = 0.2;         % COM of Pendulum
g   = 9.81;        % Gravity Constant
b   = 0.00007892;  % viscous damping at pivot of Pendulum
L   = 0.046;       % Motor Inductance
Rm  = 12.5;        % Motor Resistance
kb  = 0.031;       % Motor back emf constant
kt  = 0.031;       % Motor Torque constant
c   = 0.63;        % friction coefficient of cart


%% calculate linearization
% x_eqb = [0 pi 0 0]
AA  = I*(M+m) + M*m*(l^2);
aa  = (((m*l)^2)*g)/AA;
bb  = ((I +m*(l^2))/AA)*(c + (kb*kt)/(Rm*(r^2)));
cc  = (b*m*l)/AA;
dd  = (m*g*l*(M+m))/AA;
ee  = ((m*l)/AA)*(c + (kb*kt)/(Rm*(r^2)));
ff  = ((M+m)*b)/AA;
mm  = ((I +m*(l^2))*kt)/(AA*Rm*r);
nn  = (m*l*kt)/(AA*Rm*r);
b   = 0; % damping

% adding damping = 0.1 to the cart
A_damped    = [0.1 0 1 0; 0 0 0 1; 0 aa -bb -cc; 0 dd -ee -ff];
A_undamped  = [0 0 1 0; 0 0 0 1; 0 aa -bb -cc; 0 dd -ee -ff];
B           = [0;0; mm; nn]; 

[~,D_damped,W_damped] = eig(A_damped);
[~,D_undamped,W_undamped] = eig(A_undamped);

%% define locally stable system
A = A_damped;
k_poles           = place(A,B,[-1;-2;-3;-4]);
A_stable          = A-B*k_poles;
sys_info.A_stable = A_stable;

if(use_stable)
    [~,D,W] = eig(A_stable);
else
    % saddle?
    [~,D,W] = eig(A);
end

%% store system params
sys_info.r  = r;
sys_info.M  = M; 
sys_info.m  = m;
sys_info.I  = I; 
sys_info.l  = l; 
sys_info.g  = g;
sys_info.b  = b; 
sys_info.L  = L;
sys_info.Rm = Rm;
sys_info.kb = kb; 
sys_info.kt = kt; 
sys_info.c  = c; 

% system setup for baseline control
sys_info.A_undamped = A_undamped;
sys_info.D_undamped = D_undamped;
sys_info.W_undamped = W_undamped;
sys_info.A_damped   = A_damped;
sys_info.D_damped   = D_damped;
sys_info.W_damped   = W_damped;

% params for swing up control
sys_info.n       = 3;
sys_info.k_swing = 1.2;


% linearization params for path integral control
sys_info.A           = A;
sys_info.B           = B;
sys_info.eig_vals    = D;
sys_info.eig_vectors = W;
sys_info.x_eqb       = [0; pi; 0; 0];

% setup for local control in path integral
sys_info.use_stable     = use_stable;
sys_info.use_unstable   = use_unstable;

% local control 
sys_info.k_poles = k_poles;