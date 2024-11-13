function grad_phi = compute_gradients(x_op, phi)
    % Inputs:
    % x_op  : current operating point
    % phi   : struct of eig funs
    %
    % Outputs:
    % grad_phi  : output struct
    
    % parse dimensions
    n_dim = length(x_op);
    
    % get local grid used for phi
    xx_local    = phi.grid;
    local_axes  = phi.axis;

    % slice grid points
    axis_length = length(phi.axis{1}); % TODO: generalize to n dim
    slice = repmat({1:axis_length}, 1, n_dim);

    % Loop over each phi function to compute its gradient
    d_phi_dx = cell(1, n_dim); 
    for i = 1:n_dim
        d_phi_dx{i} = gradient(phi.phi(slice{:},i), local_axes{:});
    end  
    
    % flip direction 1 and 2 to match ndgrid with state space axis
    gradients = cell(1, n_dim);
    for i = 1:n_dim 
        order = [2, 1, 3, 4]; %TODO: generlaize to ndim
        gradients{i} = permute(d_phi_dx{i}, order);
    end
    
    % find the idx of x_op
    idx_center = num2cell(repmat(ceil(axis_length / 2), 1, n_dim));
    slice_idx = repmat({3}, 1, 4);

    % Initialize a cell array to store gradient slices
    grad_phi_x_op = nan(n_dim, n_dim);
    
    % Loop through each function and extract the gradient slice for each dimension
    for i = 1:n_dim        
        % Set the appropriate index for the current dimension (func_idx)
         slice_idx{i} = ':';
        
        % Extract the gradient slice for the current function at the specified point
        grad_phi_x_op(:,i) = gradients{i}(slice_idx{:});
    end

    % Store everything in a structured output
    grad_phi = struct();
    grad_phi.grad_phi = gradients;

    % order matters according to eigvals of A (sorted max to min)
    grad_phi.grad_phi_x_op = grad_phi_x_op;

end