function [dxdt,fx,gx] = dynamics_pendulum(x, u, sys_info)


% System parameters
m = sys_info.m;                % Mass (kg)
L = sys_info.L;                % Length of pendulum (m)
g = sys_info.g;               % Gravity (m/s^2)
c = sys_info.c;                % Damping coefficient

theta = x(1);
theta_dot = x(2);
tau = u;

theta_ddot = (tau - c * theta_dot - m * g * L * sin(theta)) / (m * L^2);
dxdt = [theta_dot; theta_ddot];

theta_ddot_no_ctrl = (- c * theta_dot - m * g * L * sin(theta)) / (m * L^2);
fx = [theta_dot; theta_ddot_no_ctrl];

theta_ddot_ctrl = 1 / (m * L^2);
gx = [0; theta_ddot_ctrl];

end