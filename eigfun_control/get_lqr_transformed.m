function lqr_params_transformed = get_lqr_transformed(D,B_transformed,Q,R)
    % inputs
    % D             - diagonal eig_vals matrix
    % B_transformed - W'*B - transformed
    % Q             - inv(W')*Q*W' - transformed
    % R             - control cost
    %
    % output
    % lqr_params (struct)
    
    % init struct
    lqr_params_transformed = struct();
    
    % calculation of LQR gain
    [K_lqr,P_lqr,e] = lqr(D,B_transformed,Q,R);
    
    % parse outputs
    lqr_params_transformed.A            = D;
    lqr_params_transformed.B            = B_transformed;
    lqr_params_transformed.Q            = Q;
    lqr_params_transformed.Q_terminal   = 10*Q;
    lqr_params_transformed.R            = R;
    lqr_params_transformed.K_lqr        = K_lqr;
    lqr_params_transformed.P_lqr        = P_lqr;
end