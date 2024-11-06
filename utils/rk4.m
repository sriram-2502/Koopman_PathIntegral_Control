function x_next = rk4(dynamics, dt_sim, x_op, u)
    x = x_op(1:2)';
    z = x_op(3:4)';
    
    k1 = z;
    states = [x;z];
    sys_info = dynamics(states,u);
    l1 = sys_info.dxdt(1:2); % take pos states only
    
    k2 = z + (dt_sim/2).*l1;
    states = [x+(dt_sim/2).*k1; z+(dt_sim/2).*l1];
    sys_info = dynamics(states,u);
    l2 = sys_info.dxdt(1:2);
    
    k3 = z + (dt_sim/2).*l2;
    states = [x+(dt_sim/2).*k2; z+(dt_sim/2).*l2];
    sys_info = dynamics(states,u);
    l3 = sys_info.dxdt(1:2);
    
    k4 = z + dt_sim.*l3;
    states = [x+dt_sim.*k3; z+dt_sim.*l3];
    sys_info = dynamics(states,u);
    l4 = sys_info.dxdt(1:2);
    
    x_update = x + (dt_sim/6).*(k1 + 2.*k2 + 2.*k3 + k4);
    z_update = z + (dt_sim/6).*(l1 + 2.*l2 + 2.*l3 + l4);
    
    x_next = [x_update;z_update];
end