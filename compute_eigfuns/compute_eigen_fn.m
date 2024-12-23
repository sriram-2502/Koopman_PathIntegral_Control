function eigen_fn = compute_eigen_fn(x_local, x_eqb, dynamics, D, W, sys_info)
    % parse inputs
    n_dim = length(x_eqb);

    % get modified linear system
    if(sys_info.use_stable)
        A_stable = sys_info.A_stable;
        dynamics_linearized = @(x,u,sys_params) A_stable*x;
    elseif(sys_info.use_unstable)
        A_unstable = sys_info.A_unstable;
        dynamics_linearized = @(x,u,sys_params) A_unstable*x;
    end

    % parse path integral setup
    if(strcmp(sys_info.id,'cart_pole'))
        path_integral_params = cart_pole_PI_params;
    elseif(strcmp(sys_info.id,'non_linear'))
        path_integral_params = nonlinear_PI_params;
    elseif(strcmp(sys_info.id,'lienar'))
        path_integral_params = linear_PI_params;
    else
        path_integral_params = PI_params;
    end
    
    % check for reverse time
    use_reverse = path_integral_params.unstable_reverse;
    if(use_reverse)
        eig_vals = -diag(D);
    else
        eig_vals = diag(D);
    end

     % make sure all eigvals are positive
     % if(any(round(diag(D))<0))
     %     disp('eigvals are negative!! cannot use compute_unstable')
     %     return
     % end

    %% open loop simualtion
    t_start = 0;
    dt_sim  = path_integral_params.dt_sim;
    t_end   = path_integral_params.t_end;
    Xout    = x_local';
    Tout    = 0;
    for t_sim = t_start:dt_sim:t_end

        % forward simulate using rk4 with no control
        x_next_full   = euler(dynamics,dt_sim,x_local,0,use_reverse,sys_info);
        x_next_linear = euler(dynamics_linearized,dt_sim,x_local,0,use_reverse,sys_info);

        % shift eqb point
        x_next_full   = x_next_full   - x_eqb;
        x_next_linear = x_next_linear - x_eqb;

        % get nonlinear part only
        x_next = x_next_full' - x_next_linear';
        
        % logs
        Tout  = [Tout;t_sim];
        Xout  = [Xout;x_next];
    end

    %% compute nonlinear part of eigfun
    integrand_convergence = cell(n_dim);
    solution_convergence  = cell(n_dim);
    for i = 1:n_dim

        % get eigval and eigvec
        lambda  = eig_vals(i);
        w       = W(:,i);

        % compute path integral
        integrand = exp(-Tout*lambda).*(w'*Xout')';
        phi_nonlinear{i} = trapz(Tout,integrand,1);
        phi_linear{i} = w'*x_local;
        phi{i} = phi_linear{i}  + phi_nonlinear{i};
    
        % check for convergence
        sol_conv = (w'*Xout')';
        nonlinear_convergence = integrand(end);
        integrand_convergence{i} = nonlinear_convergence;
        solution_convergence{i}  = sol_conv(end);
    end

    % Loop through each element in phi and assign it to phi_forward.phi
    for i = 1:n_dim
        eigen_fn.phi(i)           = phi{i};
        eigen_fn.phi_linear(i)    = phi_linear{i};
        eigen_fn.phi_nonlinear(i) = phi_nonlinear{i};
        eigen_fn.integrand(i)     = integrand_convergence{i};
        eigen_fn.sol_conv(i)      = solution_convergence{i};
    end