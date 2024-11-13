function phi = setup_path_integrals(x_op, dynamics)
% Inputs:
% x_op      : current x_op
% dynamics  : dynamics of the system
% Outputs:
% phi_local : struct
    
    %% setup grid
    [local_grid, local_axes] = setup_local_grid(x_op);
    num_elements = numel(local_grid{:});
    n_dim = length(x_op);

    % Flatten each grid into a column and concatenate into a matrix of points
    grid_points = cellfun(@(grid) grid(:), local_grid, 'UniformOutput', false);
    grid_points = [grid_points{:}];  % Concatenate into a single matrix

    % Initialize arrays to store computed values
    phi_complete  = nan(num_elements, n_dim);
    phi_linear    = nan(num_elements, n_dim);
    phi_nonlinear = nan(num_elements, n_dim); 
    phi_integrand = nan(num_elements, n_dim); 
    phi_x_op      = nan(1, n_dim);  

    %% check for stable, anti-stable or saddle
    [~, sys_info] = dynamics(x_op, 0);
    D       = sys_info.eig_vals;
    W       = sys_info.eig_vectors;
    x_eqb   = sys_info.x_eqb;

    if(all(diag(D) < 0))
        disp('---- using forward time path integrals -----')
    elseif(all(diag(D) > 0))
        disp('---- using reverse time path integrals -----')
    end

    %% Local path integral computation
    % Iterate through all grid points
    for idx = 1:num_elements
        x_local = grid_points(idx,:);

        % compute path integral at that point
        if(all(diag(D) < 0))
            phi_forward = compute_forward_time(x_local, x_eqb, dynamics, D, W);
            phi_complete(idx, 1:n_dim)   = phi_forward.phi(:)';
            phi_linear(idx, 1:n_dim)     = phi_forward.phi_linear(:)';
            phi_nonlinear(idx, 1:n_dim)  = phi_forward.phi_nonlinear(:)';
            phi_integrand(idx, 1:n_dim)  = phi_forward.integrand(:)';
            
            % Extract values for current operating point
            if(norm(x_local - x_op) <= 1e-3)
                disp('----- computing eig_fun at x_op -----')
                phi_x_op(1:n_dim) = phi_forward.phi(:)';
            end

        elseif(all(diag(D) > 0))
            % TODO: modify rk4 to solve for reverse time
            phi_reverse = compute_reverse_time(x_local, x_eqb, dynamics, D, W);
            phi_complete(idx, 1:n_dim)   = phi_reverse.phi(:)';
            phi_linear(idx, 1:n_dim)     = phi_reverse.phi_linear(:)';
            phi_nonlinear(idx, 1:n_dim)  = phi_reverse.phi_nonlinear(:)';
            phi_integrand(idx, 1:n_dim)  = phi_reverse.integrand(:)';
            
            % Extract values for current operating point
            if(norm(x_local - x_op) <= 1e-3)
                disp('----- computing eig_fun at x_op -----')
                phi_x_op(1:n_dim) = phi_reverse.phi(:)';
            end
        end
    end

    %% store everything in a struct
    phi = struct();

    % store the axis values
    phi.axis = local_axes;
    phi.grid = local_grid;
    
    % store eigfun values for each grid point
    phi.phi             = phi_complete;
    phi.phi_linear      = phi_linear;
    phi.phi_nonlinear   = phi_nonlinear;
    phi.phi_integrand   = phi_integrand;
    
    % store eigfun values at the operating point
    phi.phi_x_op           = phi_x_op;
end
