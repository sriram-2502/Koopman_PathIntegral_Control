function [dxdt, sys_info] = dynamics_linear(x, u, sys_info)
% input
% x - states
% u - control
%
% outputs
% dxdt

    
    if(sys_info.use_stable)
        dxdt    = sys_info.A_stable*x + sys_info.B*u;
    elseif(sys_info.use_unstable)
        dxdt    = sys_info.A_unstable*x + sys_info.B*u;
    else
        dxdt    = sys_info.A*x + sys_info.B*u;
    end
end






