function x_next = euler(dynamics, dt_sim, x_op, u)
    % Perform Euler integration to compute the next state
    % 
    % Inputs:
    %   dynamics   - Function handle that returns dxdt given x_op and u
    %   dt_sim     - Time step for the simulation
    %   x_op       - Current state vector (n x 1)
    %   u          - Control input (could be a vector, depending on the dynamics)
    % 
    % Output:
    %   x_next     - Next state vector (n x 1)

    % Compute the rate of change of the state
    dxdt = dynamics(x_op, u);  % Call the dynamics function to get dxdt
    
    % Update the state using Euler integration
    x_next = x_op + dt_sim * dxdt;  % Euler method: x_next = x_op + dt * dxdt
end