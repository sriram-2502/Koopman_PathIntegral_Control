function [t_riccati,P_riccati] = compute_riccati(sys_info, lqr_params_transfomed, t_span)
    
    % check controlability and observability
    A = sys_info.eig_vals;
    W = sys_info.eig_vectors;
    B = W'*sys_info.B;

    % get LQR params
    Q   = lqr_params_transfomed.Q;
    R   = lqr_params_transfomed.R;
    Q_T = lqr_params_transfomed.Q_terminal;
    
    %% Solve (numerically) for Finite-horizon LQR
    t_span_rev = flipud(t_span');
    P_initial = Q_T;
    
    options = odeset('JConstant','on', 'RelTol',1e-6, 'AbsTol',1e-6);
    [t_riccati,P_riccati]=ode45(@(t,P)riccati_ode(t,P,A,B,Q,R),t_span_rev,P_initial,options);
    P_riccati = flipud(P_riccati);
    
end