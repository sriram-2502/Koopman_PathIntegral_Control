function sys_info = linear_sys_info(x, u, sys_params)
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
    rng(102)
    A = rand(n_dim);
    B = eye(n_dim,1); % TODO: check this to make [0 0 ... 1]?
    n = length(x);
    if(~rank(ctrb(A,B))==n)
        disp('!!! system not controllable. try again')
        exit
    end
   
    % define locally stable system
    eig_vals_stable = -10*abs(rand(1,n_dim));
    K_poles_stable  = place(A,B,eig_vals_stable);
    A_stable        = A-B*K_poles_stable;
    if(~rank(ctrb(A_stable,B))==n)
        disp('!!! locally stable system not controllable. try again')
        return
    end

    % define locally unstable system
    eig_vals_unstable = 10*abs(rand(1,n_dim));
    K_poles_unstbale = place(A,B,eig_vals_unstable);
    A_unstable       = A-B*K_poles_unstbale;
    if(~rank(ctrb(A_unstable,B))==n)
        disp('!!! locally unstable system not controllable. try again')
        return
    end
 
% TODO: setup lqr option
%     Q = eye(n_dim);
%     R = eye(size(u));
%     K_lqr = lqr(A,B,Q,R);
%     A_stable = A-B*K_lqr;

% init output struct  
if(sys_params.use_stable)
    [~,D,W] = eig(A_stable);
elseif(sys_params.use_unstable)
    [~,D,W] = eig(A_unstable);
else
    % saddle?
    [~,D,W] = eig(A);
end
    
    %% parse outputs
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






