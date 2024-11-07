function u_energy = get_energy_based(dynamics, x_op)
    
    % get system info
    sys_info = dynamics(x_op, 0);

    % parse states
    x           = x_op(1);
    theta       = x_op(2);
    x_dot       = x_op(3);
    theta_dot   = x_op(4);

    % prase system params
    r   = sys_info.r;   % radius of motor shaft
    M   = sys_info.M;   % Mass of cart
    m   = sys_info.m;   % mass of pendulum
    I   = sys_info.I;   % MOI of Pendulum
    l   = sys_info.l;   % COM of Pendulum
    g   = sys_info.g;   % Gravity Constant
    b   = sys_info.b;   % viscous damping at pivot of Pendulum
    L   = sys_info.L;   % Motor Inductance
    Rm  = sys_info.Rm;  % Motor Resistance
    kb  = sys_info.kb;  % Motor back emf constant
    kt  = sys_info.kt;  % Motor Torque constant
    c   = sys_info.c;   % friction coefficient of cart

    % parse energy of the system
    E   = sys_info.E; % total energy
    Er  = sys_info.Er; % potential energy 

    % setup params for swing control
    n       = sys_info.n; 
    k_swing = sys_info.k_swing;

    % Energy based swing up control 
    accel = 2*(E-Er)*sign(theta_dot*cos(theta));
    accel = k_swing*g*(saturate_fun(accel, n*g, -n*g));

    % feedback Linearization
    u_energy = (M+m)*(accel)+ 0*x_dot-m*l*( (theta_dot)^2)*sin(theta)- m*l*(cos(theta))*( ( b*theta_dot + m*l*accel*cos(theta) + m*g*l*sin(theta) )/(I+m*l^2) );
end