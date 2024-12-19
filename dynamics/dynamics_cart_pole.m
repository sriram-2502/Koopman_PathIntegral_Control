function [dxdt, sys_info] = dynamics_cart_pole(x, u, sys_info)
% input
% x - states (pos, angle, vel, ang vel)
% u - control
%
% outputs
% sys info (struct)

% parse params
M = sys_info.M; 
m = sys_info.m;
l = sys_info.l; 
g = sys_info.g;
b = sys_info.b; 
c = sys_info.c; 
I = sys_info.I;

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

dxdt = [x_dot; theta_dot; x_ddot; theta_ddot];











