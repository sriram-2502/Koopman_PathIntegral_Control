function [t_riccati,P_riccati] = compute_riccati(sys_info, lqr_params_transfomed, t_span)
    
    % check controlability and observability
    A_transformed = sys_info.eig_vals;
    W = sys_info.eig_vectors;
    B_transformed = W'*sys_info.B;

    % get LQR params
    Q_transformed    = lqr_params_transfomed.Q;
    R_transformed    = lqr_params_transfomed.R;
    Q_T_transformed  = lqr_params_transfomed.Q_terminal;
    
    %% Solve (numerically) for Finite-horizon LQR
    t_span_rev = flipud(t_span');
    P_initial = Q_T_transformed;
    
    options = odeset('JConstant','on', 'RelTol',1e-6, 'AbsTol',1e-6);
    [t_riccati,P_riccati]=ode45(@(t,P)riccati_ode(t,P,A_transformed,B_transformed,Q_transformed,R_transformed),t_span_rev,P_initial,options);
    P_riccati = flipud(P_riccati);
    
end