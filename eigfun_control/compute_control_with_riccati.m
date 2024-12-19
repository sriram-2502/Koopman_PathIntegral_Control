function u = compute_control_with_riccati(lqr_params_transformed, sys_info, phi_x_op, grad_phi_x_op, t_span_curr, x_op)
    
    % check if system is linear or not
    if(isempty(grad_phi_x_op))
        use_linear = true;
    else
        use_linear = false;
    end

    % parse stuff
    B = sys_info.dynamics_g(x_op(1),x_op(2),x_op(3),x_op(4));
    D = sys_info.eig_vals;
    Q = lqr_params_transformed.Q;
    R = lqr_params_transformed.R;
    W = grad_phi_x_op; % linearized eig fun at x_op
    
    A_transformed           = D;
    B_transformed           = W'*B;
    Q_transformed           = inv(W)*Q*inv(W');
    lqr_params_transformed  = get_lqr(A_transformed,B_transformed,Q_transformed,R);

   % solve riccati at curren time step
   [t_riccati,P_riccati] = compute_riccati(lqr_params_transformed, t_span_curr);
   P_riccati_curr = reshape(P_riccati(1,:),size(lqr_params_transformed.A));

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