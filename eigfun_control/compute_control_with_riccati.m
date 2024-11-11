function u = compute_control_with_riccati(x_op, sys_info, lqr_params, phi, grad_phi, t_span_curr)

    % get riccati solution for current time step 
    A = sys_info.A_koopman;
    W = sys_info.W;
    g = W'*sys_info.g(x_op'); % eval g(x) at x_op
    Q = inv(W')*lqr_params.Q*W';
    R = lqr_params.R;

    [t_riccati,P_riccati] = compute_riccati(sys_info, lqr_params, t_span_curr);
    P_riccati_curr = reshape(P_riccati(1,:),size(A));

    % get eig fun and its grad at current x_op
    phi = phi.phi_x_op;
    grad_phi = grad_phi.grad_phi_x_op;

    % get control
    u = -inv(R)*g'*(phi'*P_riccati_curr*grad_phi)';
    
end