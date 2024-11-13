function [dxdt, sys_info] = dynamics_linear(x, u)
% input
% x - states
% u - control
%
% outputs
% sys info (struct)

% init output struct
sys_info = struct();
n_dim = length(x);

    %% Compute dxdt
    A = rand(n_dim);
    B = eye(n_dim,1);
    
    if(~rank(ctrb(A,B))==rank(A))
        disp('!!! system not controllable. try again')
        exit
    end
    [~,D,W] = eig(A);
    
    dxdt = A*x + B*u;
    
    sys_info.dxdt        = dxdt;
    sys_info.A           = A;
    sys_info.B           = B;
    sys_info.eig_vals    = D;
    sys_info.eig_vectors = W;
    sys_info.x_eqb       = zeros(size(x));
   
    %% define locally stable system
    K_poles = place(A,B,-(1:n_dim));
    A_stable = A-B*K_poles;
    sys_info.A_stable = A_stable;
    
    [~,D,W] = eig(A_stable);
    sys_info.eig_vals    = D;
    sys_info.eig_vectors = W;

end






