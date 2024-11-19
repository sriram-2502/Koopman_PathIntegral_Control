function [dxdt, sys_info] = dynamics_linear(x, u, sys_params)
% input
% x - states
% u - control
%
% outputs
% sys info (struct)

% init output struct
sys_info = struct();
n_dim = length(x);

if(nargin<=2)
    use_stable   = true;
    use_unstable = false;
else
    use_stable   = sys_params.use_stable;
    use_unstable = sys_params.use_unstable;
end

    %% Compute dxdt
    rng(200)
    A = rand(n_dim);
    B = eye(n_dim,1); % TODO: check this to make [0 0 ... 1]
    n = length(x);
    
    if(~rank(ctrb(A,B))==n)
        disp('!!! system not controllable. try again')
        exit
    end
   
    % define locally stable system
    K_poles_stable  = place(A,B,1:n_dim);
    A_stable        = A-B*K_poles_stable;

    % define locally unstable system
    K_poles_unstbale = place(A,B,1:n_dim);
    A_unstable       = A-B*K_poles_unstbale;
 
% TODO: setup lqr option
%     Q = eye(n_dim);
%     R = eye(size(u));
%     K_lqr = lqr(A,B,Q,R);
%     A_stable = A-B*K_lqr;
    
    if(use_stable)
        [~,D,W] = eig(A_stable);
        dxdt    = A_stable*x + B*u;
    elseif(use_unstable)
        [~,D,W] = eig(A_unstable);
        dxdt    = A_unstable*x + B*u;
    else
        % saddle?
        [~,D,W] = eig(A);
        dxdt    = A*x + B*u;
    end

    
    %% parse outputs
    sys_info.dxdt           = dxdt;
    sys_info.A              = A;
    sys_info.A_stable       = A_stable;
    sys_info.A_unstable     = A_unstable;
    sys_info.B              = B;
    sys_info.eig_vals       = D;
    sys_info.eig_vectors    = W;
    sys_info.x_eqb          = zeros(size(x));
    sys_info.use_stable     = use_stable;
    sys_info.use_unstable   = use_unstable;

end






