function phi_reverse = compute_reverse_time(x_local, x_eqb, dynamics, D, W)
    
    % parse inputs
    n_dim = length(x_eqb);
    sys_info = dynamics(x_local,0);
    A = sys_info.A;
    dynamics_linear = @(x,u)A*x;

    %% open loop simualtion
    t_start = 0;
    dt_sim  = 0.01;
    t_end   = 10;
    Xout    = x_local;
    for t_sim = t_start:dt_sim:t_end

        % forward simulate using rk4 with no control
        x_next_full   = rk4(dynamics,dt_sim,x_local,0);
        x_next_linear = rk4(dynamics_linear,dt_sim,x_local,0);

        % shift eqb point
        x_next_full   = x_eqb - x_next_full;
        x_next_linear = x_eqb - x_next_linear;

        % get nonlinear part only
        x_next = x_next_full - x_next_linear;
        
        % logs
        Tout  = [Tout;t_sim];
        Xout  = [Xout; x_next];
    end

    %% compute nonlinear part of eigfun
    eig_vals = diag(D);
    for i = 1:n_dim

        % get eigval and eigvec
        lambda  = eig_vals(i);
        w       = W(:,i);

        % compute path integral
        integrand = exp(Tout(:)*lambda).*Xout;
        phi_nonlinear{i} = trapz(Tout,integrand,1);
        phi_linear{i} = w'*x_i;
        phi{i} = phi_linear{i}  + phi_nonlinear{i};
    
        % check for convergence (use abs value)
        abs_integrand = exp(Tout(:)*lambda).*abs(Xout);
        integrand{i} = abs_integrand(end);
    end

    phi_reverse.phi{1:n_dim} = phi{:};
    phi_reverse.phi_linear{1:n_dim} = phi_nonlinear{:};
    phi_reverse.phi_linear{1:n_dim} = phi_nonlinear{:};
    phi_reverse.integrand{1:n_dim} = integrand{:};