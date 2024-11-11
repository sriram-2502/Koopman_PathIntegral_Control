function phi = setup_path_integrals(x_op, dynamics)
% Inputs:
% x_op      : current x_op
% dynamics  : dynamics of the system
% Outputs:
% phi_local : struct 
    
    %% setup grid
    n_dim = length(x_op);
    [xx_local, local_axes] = setup_local_grid(x_op);
    % Unpack the grid
    axis1 = xx_local{1}; %transpose for consistency?
    axis2 = xx_local{2};
    axis3 = xx_local{3};
    axis4 = xx_local{4};

    %% check for stable, anti-stable or saddle
    % TODO: add local controllers to make system locally stable
    [~,sys_info] = dynamics(x_op,0);
    D       = sys_info.eig_vals;
    W       = sys_info.eig_vectors;
    x_eqb   = sys_info.x_eqb;

    use_forward_time    = [];
    use_reverse_time    = [];
    use_forward_reverse = [];

    if(all(diag(D)<0))
        use_forward_time = true;
    elseif(all(diag(D)>0))
        use_reverse_time = true;
    else
        use_forward_reverse = true;
    end

    %% Local path integral computation
    for i = 1:size(axis1,1)
        i
        for j = 1:size(axis2,2)
            for k = 1:size(axis3,3)
                for l = 1:size(axis4,4)

                    % get local grind point
                    x_local = [axis1(i,j,k,l); axis2(i,j,k,l); axis3(i,j,k,l); axis4(i,j,k,l)];
                    
                    % compute path inegral at that point
                    if(use_forward_time)
                        phi = compute_forward_time(x_local, x_eqb, dynamics, D, W);
        
                    elseif(use_reverse_time)
                        % TODO: modify rk4 to solve for reverse time
                        phi = compute_reverse_time(x_local, x_eqb, dynamics, D, W);
                    end

                end
            end
        end
    end
    
    %% Extract path integrals at x_op

    % Calculate the squared distance between each grid point and x_op
    distances = (axis1 - x_op(1)).^2 + (axis2 - x_op(2)).^2 + (axis3 - x_op(3)).^2 + (axis4 - x_op(4)).^2;
    
    % Find the index of the grid point with the minimum distance
    [~, idx_min] = min(distances(:));
    
    % Convert the linear index to subscripts (3D indices)
    [idx_x1, idx_x2, idx_x3, idx_x4] = ind2sub(size(distances), idx_min);
    phi_x_op = cell(1,n_dim);
    for i = 1:n_dim
        phi_x_op{i} = phi.phi{i}(idx_x1, idx_x2, idx_x3, idx_x4);
    end

    %% store everything in a struct
    phi = struct();
    phi.local_grid.axis{1:n_dim} = xx_local{:};
    phi.phi{1:n_dim} = phi.phi{:};

    % save the integrand in the same struct
    phi.phi_integrand{1:n_dim} = phi.integrand{:};

    % order matters according to eigvals of A (sorted max to min)
    phi.phi_x_op = phi_x_op{:};
end