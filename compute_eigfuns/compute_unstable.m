function phi_forward = compute_unstable(x_local, x_eqb, dynamics, D, W, sys_info)
    
    % parse inputs
    n_dim               = length(x_eqb);
    A_unstable          = sys_info.A_unstable;
    dynamics_linearized = @(x,u,sys_params) A_unstable*x;
    
    % check for reverse time
    use_reverse = true;

     % make sure all eigvals are positive
     % if(any(round(diag(D))<0))
     %     disp('eigvals are negative!! cannot use compute_unstable')
     %     return
     % end

    %% open loop simualtion
    t_start = 0;
    dt_sim  = 0.01;
    t_end   = 2;
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
    if(use_reverse)
         eig_vals = -diag(D);
    else
        eig_vals = diag(D);
    end
    integrand_nonlinear = cell(n_dim);
    for i = 1:n_dim

        % get eigval and eigvec
        lambda  = eig_vals(i);
        w       = W(:,i);

        % compute path integral
        integrand = exp(-Tout*lambda).*(w'*Xout')';
        phi_nonlinear{i} = trapz(Tout,integrand,1);
        phi_linear{i} = w'*x_local;
        phi{i} = phi_linear{i}  + phi_nonlinear{i};
    
        % check for convergence (use abs value)
        integrand_nonlinear{i} = integrand(end);
    end

    % Loop through each element in phi and assign it to phi_forward.phi
    for i = 1:n_dim
        phi_forward.phi(i)           = phi{i};
        phi_forward.phi_linear(i)    = phi_linear{i};
        phi_forward.phi_nonlinear(i) = phi_nonlinear{i};
        phi_forward.integrand(i)     = integrand_nonlinear{i};
    end