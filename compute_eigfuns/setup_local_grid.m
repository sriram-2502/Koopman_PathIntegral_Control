function [xx_local, local_axes] = setup_local_grid(x_op)
    % setup_local_grid_nd creates an n-dimensional local grid centered around x_op.
    %
    % Inputs:
    %   x_op        - Center point of the grid (n x 1 vector)
    %
    % Outputs:
    %   xx_local    - Cell array containing meshgrid values for each dimension
    %   local_axes  - Cell array of grid points along each axis

    % grid params
    grid_size = 0.02;   % Size of neighborhood around x_op
    step_size = 0.01;   % Step size for the local grid
    
    n_dim = length(x_op);  % Number of dimensions

    % Initialize local_axes to store grid points for each dimension
    local_axes = cell(n_dim, 1);
    
    % Create grid points for each dimension
    for i = 1:n_dim
        local_axes{i} = x_op(i) - grid_size : step_size : x_op(i) + grid_size;
    end
    
    % Generate meshgrid based on local_axes for n-dimensional grid
    [xx_local{1:n_dim}] = ndgrid(local_axes{:});
end