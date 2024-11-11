function u = compute_control(sys_info, lqr_params, P_riccati, phi, grad_phi)

    % parse system info
    W = sys_info.W;
    B = W'*sys_info.B; %lienarized version of g(x)

    % parse lqr params
    R = lqr_params.R;

    % get eig fun and its grad at current x_op
    phi = phi.phi_x_op;
    grad_phi = grad_phi.grad_phi_x_op;

    % get control
    u = -inv(R)*B'*(phi'*P_riccati*grad_phi)';
    
end