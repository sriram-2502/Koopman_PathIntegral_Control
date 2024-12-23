function dxdt = dynamics_saddle(x,u, sys_info)

dxdt = sys_info.f(x) + sys_info.g(x)*u;
