function u = compute_control(sys_info, lqr_params, phi_x_op, grad_phi_x_op)

    % parse system info
    B           = sys_info.B; %lienarized version of g(x)
    P_riccati   = lqr_params.P_riccati_curr;

    % parse lqr params
    R = lqr_params.R;

    % get control
    u = -inv(R)*B'*(phi_x_op*P_riccati*grad_phi_x_op)';
    
end