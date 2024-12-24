function dxdt = dynamics_nonlinear(x,u, sys_info)

dxdt = sys_info.dynamics_f(x) + sys_info.dynamics_g(x)*u;
