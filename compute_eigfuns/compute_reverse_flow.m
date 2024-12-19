function phi_reverse = compute_reverse_flow(x_local, x_eqb, dynamics, D, W, sys_info)
    
    % parse inputs
    n_dim           = length(x_eqb);
    A_unstable      = sys_info.A_unstable;
    dynamics_linear = @(x,u,sys_params)A_unstable*x;

     % check for reverse time
    if(sys_info.use_unstable)
        use_reverse = true; 
    elseif(sys_info.use_stable)
        use_reverse = false;
    end

     % make sure all eigvals are negative
     if(any(ceil(diag(D))<0))
         disp('eigvals are positive!! cannot use reverse_time')
         return
     end

    %% open loop simualtion
    t_start = 0;
    dt_sim  = 0.01;
    t_end   = 1;
    Xout    = x_local';
    Tout    = 0;
    for t_sim = t_start:dt_sim:t_end

        % forward simulate using rk4 with no control
        x_next_full   = euler(dynamics,dt_sim,x_local,0,use_reverse, sys_info);
        x_next_linear = euler(dynamics_linear,dt_sim,x_local,0,use_reverse,sys_info);

        % shift eqb point
        x_next_full   = x_next_full   - x_eqb;
        x_next_linear = x_next_linear - x_eqb;

        % get nonlinear part only
        x_next = x_next_full' - x_next_linear';
        
        % logs
        Tout  = [Tout;t_sim];
        Xout  = [Xout; x_next];
    end

    %% compute nonlinear part of eigfun
    eig_vals = diag(D);
    integrand_nonlinear = cell(n_dim);
    for i = 1:n_dim

        % get eigval and eigvec
        lambda  = eig_vals(i);
        w       = W(:,i);

        % compute path integral
        integrand = exp(Tout*lambda).*(w'*Xout')';
        phi_nonlinear{i} = trapz(Tout,integrand,1);
        phi_linear{i} = w'*x_local;
        phi{i} = phi_linear{i}  + phi_nonlinear{i};
    
        % check for convergence (use abs value)
        abs_integrand = exp(Tout*lambda).*abs(Xout);
        integrand_nonlinear{i} = abs_integrand(end);
    end

    % Loop through each element in phi and assign it to phi_forward.phi
    for i = 1:n_dim
        phi_reverse.phi(i)           = phi{i};
        phi_reverse.phi_linear(i)    = phi_linear{i};
        phi_reverse.phi_nonlinear(i) = phi_nonlinear{i};
        phi_reverse.integrand(i)     = integrand_nonlinear{i};
    end