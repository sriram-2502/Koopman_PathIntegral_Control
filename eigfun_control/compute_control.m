function u = compute_control(lqr_params_transformed, P_riccati_curr, phi_x_op, grad_phi_x_op)

    % parse system info
    B = lqr_params_transformed.B; %lienarized version of g(x) in transformed coordinates

    % parse lqr params
    R = lqr_params_transformed.R;

    % get control
    u = -inv(R)*B'*(phi_x_op*P_riccati_curr*grad_phi_x_op)';
    
end