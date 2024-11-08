function phi_local = setup_path_integrals(x_op, dynamics)
% Inputs:
% x_op      : current x_op
% dynamics  : dynamics of the system
% Outputs:
% phi_local : struct 
    
    %% setup grid
    [xx_local, local_axes] = setup_local_grid(x_op);
    % Unpack the grid
    axis1 = xx_local{1}; %transpose for consistency?
    axis2 = xx_local{2};
    axis3 = xx_local{3};
    axis4 = xx_local{4};

    %% Init eigfuns
    phi1 = nan(size(axis1)); 
    phi2 = nan(size(axis2));
    phi3 = nan(size(axis3));
    phi4 = nan(size(axis4));
    nl_phi1 = nan(size(phi1)); 
    nl_phi2 = nan(size(phi2));
    nl_phi3 = nan(size(phi3));
    nl_phi4 = nan(size(phi4));
    nl_phi1_integrand = nan(size(nl_phi1)); 
    nl_phi2_integrand = nan(size(nl_phi2)); 
    nl_phi3_integrand = nan(size(nl_phi3)); 
    nl_phi4_integrand = nan(size(nl_phi4)); 

    %% check for stable, anti-stable or saddle
    % TODO: add local controllers to make system locally stable
    sys_info = dynamics(x_op,0);
    D = sys_info.eig_vals;
    W = sys_info.eig_vectors;
    if(all(diag(D)>0))
        use_forward_time = true;
    elseif(all(diag(D)<0))
        use_reverse_time = true;
    else
        use_forward_reverse = true;
    end

    %% Local path integral computation
    for i = 1:size(axis1,1)
        for j = 1:size(axis2,2)
            for k = 1:size(axis3,3)

                % get local grind point
                x_local = [axis1(i,j,k); axis2(i,j,k); axis3(i,j,k)];
                
                % compute path inegral at that point
                if(use_forward_time)
                    phi_forward = compute_forward_time(x_local, x_eqb, dynamics, D, W);
                    phi1(i,j,k) = phi_forward.phi.phi;
                    nl_phi1(i,j,k) = phi_forward.nl_phi;
                    nl_phi1_integrand(i,j,k) = phi_forward.integrand;
    
                elseif(use_reverse_time)
                    % TODO: modify rk4 to solve for reverse time
                    phi_reverse = compute_reverse_time(x_local, x_eqb, dynamics, D, W);
                    phi1(i,j,k) = phi_reverse.phi.phi;
                    nl_phi1(i,j,k) = phi_reverse.nl_phi;
                    nl_phi1_integrand(i,j,k) = phi_reverse.integrand;
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
    phi1_x_op = phi1(idx_x1, idx_x2, idx_x3, idx_x4);
    phi2_x_op = phi2(idx_x1, idx_x2, idx_x3, idx_x4);
    phi3_x_op = phi3(idx_x1, idx_x2, idx_x3, idx_x4);
    phi4_x_op = phi4(idx_x1, idx_x2, idx_x3, idx_x4);

    %% store everything in a struct
    phi_local = struct();
    phi_local.local_grid.xx1_local = axis1;
    phi_local.local_grid.xx2_local = axis2;
    phi_local.local_grid.xx3_local = axis3;
    phi_local.local_grid.xx4_local = axis4;
    phi_local.phi1 = phi1;
    phi_local.phi2 = phi2;
    phi_local.phi3 = phi3;
    phi_local.phi4 = phi4;

    % save phi at the current operating point
    phi_local.phi1_x_op = phi1_x_op;
    phi_local.phi2_x_op = phi2_x_op;
    phi_local.phi3_x_op = phi3_x_op;
    phi_local.phi4_x_op = phi4_x_op;
    
    % save the integrand in the same struct
    phi_local.phi1_integrand = nl_phi1_integrand;
    phi_local.phi2_integrand = nl_phi2_integrand;
    phi_local.phi3_integrand = nl_phi3_integrand;
    phi_local.phi4_integrand = nl_phi4_integrand;

    % order matters according to eigvals of A (sorted max to min)
    phi_local.phi_x_op = [phi1_x_op; phi2_x_op; phi3_x_op; phi4_x_op];
end