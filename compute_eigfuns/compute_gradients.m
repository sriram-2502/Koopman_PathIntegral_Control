function grad_phi = compute_gradients(x_op, phi)
    % Inputs:
    % x_op  : current operating point
    % phi   : struct of eig funs
    %
    % Outputs:
    % grad_phi  : output struct

    n = numel(phi.phi);  % Get the number of elements in phi.phi
    % Loop to assign each element to dynamically named variables
    for i = 1:n
        eval(sprintf('phi%d = phi.phi{%d};', i, i));
    end

    [xx_local, local_axes] = setup_local_grid_nd(x_op);
    dims = numel(xx_local);  % Number of dimensions
    % Unpack the grid
    for i = 1:n
        eval(sprintf('xx%d_local =xx_local{%d};', i, i));
        eval(sprintf('local_x%d = local_axes{%d};', i, i));
    end

    % Initialize cell arrays to store gradients
    d_phi_dx = cell(n, dims); 
    % Loop over each phi function to compute its gradient
    for i = 1:n
        % Compute gradient for phi{i}
        [gradients{1:dims}] = gradient(phi.phi{i}, local_x{:});
        
        % Store each gradient component in the corresponding cell
        for j = 1:dims
            d_phi_dx{i, j} = gradients{j};
        end
    end

    % Initialize a cell array to store gradient vectors for each phi
    grad_phi = cell(1, n);
    
    % Loop through each phi and store the gradients in a specified order
    for i = 1:n
        % Define the order of dimensions (e.g., [2, 1, 3] for a custom order)
        order = [2, 1, 3, 4];  % Adjust the order as needed for n dimensions
        
        % Store the gradients in the specified order
        grad_phi{i} = cell(1, dims);
        for j = 1:dims
            grad_phi{i}{j} = d_phi_dx{i, order(j)};
        end
    end

    % Sum the squared differences across all dimensions
    for dim = 1:n
        distances = distances + (xx_local{dim} - x_op(dim)).^2;
    end
    
    %Find the index of the grid point with the minimum distance
    [~, idx_min] = min(distances(:));
    
    % Convert the linear index to subscripts (3D indices)
    % Initialize a cell array to hold the indices
    indices = cell(1, n);
    [indices{:}] = ind2sub(size(distances), idx_min);
    
    % Convert cell array to numeric array for easier handling
    indices_num = cell2mat(indices);  % 1 x n vector

    % Initialize a cell array to store gradient vectors for each phi at closest point
    grad_phi_x_op = cell(1, n);

    % Extract gradients at the point closest to x_op
    for i = 1:n
        grad_phi_x_op{i} = zeros(1, dims);
        for j = 1:dims
            grad_phi_x_op{i}(j) = d_phi_dx{i, j}(indices_num{:});
        end
    end

    % Store everything in a structured output
    grad_phi = struct();
    grad_phi.grad_phi{1:n_dim} = grad_phi{:}{:};

    % order matters according to eigvals of A (sorted max to min)
    grad_phi.grad_phi_x_op = grad_phi_x_op{:}{:};

end