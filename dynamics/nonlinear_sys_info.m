function sys_info = nonlinear_sys_info(sys_params)
%% system description
% nonlinear ode x_dot = f(x) + Bu
n = 2; m = 1;
x = sym('x',[n,1],'real');

% define analytical eig_funs
Lam = [2;1];
phiST_ana = @(x1,x2)x1-2*x2;
phiUS_ana = @(x1,x2)x1+sin(x2);

Phi_ana = @(x1,x2)[phiUS_ana(x1,x2); phiST_ana(x1,x2)];
phi_x = Phi_ana(x(1),x(2));
dPhi_dx = simplify(jacobian(phi_x,x));

% get system dynamics
f_x = inv(dPhi_dx)*diag(Lam)*phi_x;

if(sys_params.use_linear_riccati)
    g_x = [cos(x(2)); 0]; % test for linearized riccati case
else
    g_x = simplify(inv(dPhi_dx)*[1;2]);
end

eig_scale = 0.001;
f_sim = inv(dPhi_dx)*diag(eig_scale*Lam)*phi_x;
f_pi = matlabFunction(f_sim,'vars',{x}); 
f = matlabFunction(f_x,'vars',{x}); 
g = matlabFunction(g_x,'vars',{x}); 
f_control = @(t,x,u) f(x) + g(x)*u; 

% grad phi of unstable nonlinear eig fun
grad_phiUS_ana = gradient(phiUS_ana(x(1),x(2)),x);
dphiUS_ana_x1 = @(x1,x2) ones(size(x1));
dphiUS_ana_x2 = @(x1,x2) cos(x2);

grad_phisT_ana = gradient(phiST_ana(x(1),x(2)),x);
dphiST_ana_x1 = @(x1,x2) ones(size(x1)); 
dphiST_ana_x2 = @(x1,x2) -2.*ones(size(x1)); 


%% Equilibrium point/s, linearization and non-linear part
xEq_struct = solve(f_x==0);
xEq = [xEq_struct.x1(1);xEq_struct.x2(1)];
A = eval(subs(jacobian(f_x),[x(1) x(2)]',xEq));
B = eval(subs((g_x),[x(1) x(2)]',xEq));
[~,D,W] = eig(A);
[dVal,dIdx] = sort(diag(D),'descend');
% arrange D and W in {unstable|stable} cofig
D = diag([D(dIdx(1),dIdx(1)),D(dIdx(2),dIdx(2))]);
W = [W(:,dIdx(1)),W(:,dIdx(2))];
% un-stable part
evUS = D(1,1);
wUS = W(:,1);
% stable part
evST = D(2,2);
wST = W(:,2);

% scale w1 and w2 to match the linearization of the analytical
% eigenfunction values
wUS = wUS./min(wUS(wUS~=0));
wST = wST./min(wST(wST~=0));
W = [wUS,wST];
% define nonlinear part x_dot = Ax + fn(x)
fn = f_x - A*x;

% define matlab functions
wUSfn = matlabFunction(wUS'*fn,'vars',{x(1), x(2)});
grad_fn_US = wUS'*simplify(jacobian(fn,x));
wUSfn_grad_x1 = matlabFunction(grad_fn_US(1),'vars',{x(1), x(2)});
wUSfn_grad_x2 = matlabFunction(grad_fn_US(2),'vars',{x(1), x(2)});

wSTfn = matlabFunction(wST'*fn,'vars',{x(1), x(2)});
grad_fn_ST = wST'*simplify(jacobian(fn,x)');
wSTfn_x1 = matlabFunction(grad_fn_ST(1),'vars',{x(1), x(2)});
wSTfn_x2 = matlabFunction(grad_fn_ST(2),'vars',{x(1), x(2)});

% get analytical controller
phi_analytical_function = @(x)[phiUS_ana(x(1),x(2)); phiST_ana(x(1),x(2))];
grad_phi_analytical_function = matlabFunction(dPhi_dx,'vars',{x});

%% setup system parameters for xdot = f(x) + B*u
sys_info.A              = A;
sys_info.B              = B;
sys_info.A_koopman      = D;
sys_info.f              = f;
sys_info.g              = g;
sys_info.state_dim      = n;
sys_info.ctrl_dim       = m;
sys_info.eig_vectors    = W;
sys_info.eig_vals       = D;
sys_info.x_eqb          = zeros(n,1);
sys_info.A_unstable     = A;
sys_info.A_stable       = A;
sys_info.id             = "non_linear";
sys_info.eigen_fun      = Phi_ana;
