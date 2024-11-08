function sys_info = dynamics_cart_pole(x, u)
% input
% x - states (pos, angle, vel, ang vel)
% u - control
%
% outputs
% sys info (struct)

% init output struct
sys_info = struct();

% parameters
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

% parse states
theta     = x(2);
x_dot     = x(3);
theta_dot = x(4);

% parse control
F = u;

%% Compute dxdt

alpha_a = ((m^2)*(l^2)*((sin(theta))^2)+ M*m*l^2 +(M+m)*I);
x_ddot  = (b*m*l*theta_dot*cos(theta) + (m^2)*(l^2)*g*sin(theta)*cos(theta) + (I + m*(l^2))*(F-c*x_dot+ m*l*sin(theta)*theta_dot^2) )/alpha_a;
theta_ddot = -(F*m*l*cos(theta)-c*m*l*x_dot*cos(theta) + (m^2)*(l^2)*(theta_dot^2)*sin(theta)*cos(theta)+ (M+m)*(b*theta_dot + m*g*l*sin(theta)))/alpha_a;

dxdt = [x_ddot; theta_ddot];
sys_info.dxdt = dxdt;

%% calculate linearization

AA  = I*(M+m) + M*m*(l^2);
aa  = (((m*l)^2)*g)/AA;
bb  = ((I +m*(l^2))/AA)*(c + (kb*kt)/(Rm*(r^2)));
cc  = (b*m*l)/AA;
dd  = (m*g*l*(M+m))/AA;
ee  = ((m*l)/AA)*(c + (kb*kt)/(Rm*(r^2)));
ff  = ((M+m)*b)/AA;
mm  = ((I +m*(l^2))*kt)/(AA*Rm*r);
nn  = (m*l*kt)/(AA*Rm*r);
A   = [0 0 1 0; 0 0 0 1; 0 aa -bb -cc; 0 dd -ee -ff];
B   = [0;0; mm; nn]; 

sys_info.A = A;
sys_info.B = B;

[~,D,W] = eig(A);
sys_info.eig_vals    = D;
sys_info.eig_vectors = W;

% TODO: check if A is the same even if eqb point is [0 pi 0 0]

%% setup params for energy shaping
% Total energy
E = m*g*l*(1-cos(theta)) + (1/2)*(I + m*l^2)*(theta_dot^2);

% Potential Energy
Er  = 2*m*g*l;

sys_info.E = E;
sys_info.Er = Er;

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

% params for swing up control
sys_info.n       = 3;
sys_info.k_swing = 1.2;






