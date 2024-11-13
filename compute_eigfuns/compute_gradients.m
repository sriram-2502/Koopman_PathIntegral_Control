function grad_phi_x_op = compute_gradients(x_op, phi)
% Compute the gradient of phi at the operating point x_op.
%
% Inputs:
% x_op     : current operating point (1 x n_dim array)
% phi      : struct containing phi values, phi_linear, phi_nonlinear, and phi_integrand
%
% Outputs:
% grad_phi_x_op : Matrix containing the gradient of phi at x_op for each dimension

    n_dim = length(x_op);   % Dimension of the space
    local_grid = phi.grid;  % Local grid values
    local_axes = phi.axis;  % Axis values for each dimension

    % Flatten each grid into a column and concatenate into a matrix of points
    grid_points = cellfun(@(grid) grid(:), local_grid, 'UniformOutput', false);
    grid_points = [grid_points{:}];  % Concatenate into a single matrix
    
    % Initialize variables to store each dimension's phi values
    grad_phi_x_op = cell(n_dim, n_dim);
    
    % loop for each eig fun
    for dim = 1:n_dim
        % Get the grid for the current dimension (axis values for this dimension)
        grid  = local_axes{dim};
        delta = grid(2) - grid(1);  % Step size (assumed uniform)
        
        % Find the index corresponding to the center (midpoint) of the grid
        idx_center = ceil(length(grid) / 2);  % Center index of the grid

        % Extract the phi values for this dimension (phi_dim)
        phi_dim = phi.phi{dim};  % phi values for the current dimension
        
        % Initialize the gradient vector at the center
        grad_phi = nan(n_dim, 1);  % Initialize gradient vector for the current dimension
        
        % loop for each direction
        for dir = 1:length(x_op)
            % Compute the central difference formula in the current direction
            % TODO: generalize to n dim
            grad_phi(dir) = (phi_dim(idx_center + 1, dir) - phi_dim(idx_center - 1, dir)) / (2 * delta);
        end
        
        % Store the gradient for the current dimension
        grad_phi_x_op{dim} = grad_phi;  % Store the gradient vector in the cell array
    end
    
end
