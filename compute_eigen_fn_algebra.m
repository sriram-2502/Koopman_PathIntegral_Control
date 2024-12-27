function eigen_fn_algebra = compute_eigen_fn_algebra(x_local, x_eqb, dynamics, D,D_negative, W, sys_info, k)
    % parse inputs
    n_dim = length(x_eqb);

    % get modified linear system
    if(sys_info.use_stable)
        A_stable = sys_info.A_stable;
        dynamics_linearized = @(x,u,sys_params) A_stable*x;
    elseif(sys_info.use_unstable)
        A_unstable = sys_info.A_unstable;
        dynamics_linearized = @(x,u,sys_params) A_unstable*x;
    else
        A = sys_info.A;
        dynamics_linearized = @(x,u,sys_params) A*x;
    end

    % parse path integral setup
    if(strcmp(sys_info.id,'cart_pole'))
        path_integral_params = cart_pole_PI_params;
    elseif(strcmp(sys_info.id,'non_linear'))
        path_integral_params = nonlinear_PI_params;
    elseif(strcmp(sys_info.id,'linear'))
        path_integral_params = linear_PI_params;
    elseif(strcmp(sys_info.id,'pendulum'))
        path_integral_params = pendulum_PI_params;
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

     % setup for algebra
     idx1            = ismember(diag(D), D_negative);
     [lambda2, idx2] = max(diag(D)); % find max positive eig val
     lambda1         = D_negative; % negative eig val
     w1              = W(:,idx1);
     w2              = W(:,idx2);
     Q               = w1*w2';
     lambda          = lambda1 + k*lambda2; % make sure lambda is always positive

    %% open loop simualtion
    t_start = 0;
    dt_sim  = path_integral_params.dt_sim;
    t_end   = path_integral_params.t_end;
    Xout    = x_local';
    Tout    = 0;
    Xout_linear = x_local';
    Xout_full   = x_local';
    
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
        Xout_full   = [Xout_full; x_next_full'];
        Xout_linear = [Xout_linear; x_next_linear'];
    end

    %% compute nonlinear part of eigfun

    % compute path integral
    integrand = exp(-Tout*lambda).*((w2'*Xout_full').^(k-1)*(Xout_full)*(Q+k*Q')*(Xout'))';
    phi_nonlinear = trapz(Tout,integrand,1);
    phi_linear = (w1'*x_local)*(w2'*x_local)^k;
    phi = phi_linear + phi_nonlinear;

    % check for convergence
    sol_conv = (w2'*Xout_full').^(k-1)*(Xout_full)*(Q+k*Q')*(Xout');
    integrand_convergence = integrand(end);
    solution_convergence  = sol_conv(end);


    % Loop through each element in phi and assign it to phi_forward.phi
    eigen_fn_algebra.phi           = phi;
    eigen_fn_algebra.phi_linear    = phi_linear;
    eigen_fn_algebra.phi_nonlinear = phi_nonlinear;
    eigen_fn_algebra.integrand     = integrand_convergence;
    eigen_fn_algebra.sol_conv      = solution_convergence;
