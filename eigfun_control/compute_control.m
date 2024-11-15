function u = compute_control(lqr_params_transformed, P_riccati_curr, phi_x_op, grad_phi_x_op)
    
    % check if system is linear or not
    if(isempty(grad_phi_x_op))
        use_linear = true;
    else
        use_linear = false;
    end

    % parse system info
    B = lqr_params_transformed.B; %lienarized version of g(x) in transformed coordinates

    % parse lqr params
    R = lqr_params_transformed.R;

    % get control
    if(use_linear)
        % if linear, no need for grad_phi
        u = -inv(R)*B'*(P_riccati_curr*phi_x_op');
    else
        u = -inv(R)*B'*(phi_x_op*P_riccati_curr*grad_phi_x_op)';
    end

    
end