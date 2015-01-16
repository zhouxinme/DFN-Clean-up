%% Crank-Nicolson Eqns for Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura
% 
% Update Log by zhouxinme:
% Jan 15th 2014: Clean up the DFN model.

function [x_nxtf, z_nxtf, varargout] = cn_dfn_zhouxinme2(x,z,Cur_vec,p)

Cur_prv = Cur_vec(1);
Cur = Cur_vec(2);
Cur_nxt = Cur_vec(3);

%% Parameters

% Newton params
maxIters = 100;
tol = 1e-6;

% Sparse Linear System Solver Params
linEqnTol = 1e-6; % Default is 1e-6
linEqnMaxIter = 100;  % Default is 20

% Length of states x and z
Nx = length(x);
Nz = length(z);

% Preallocate
x_nxt = zeros(Nx,maxIters);
z_nxt = zeros(Nz,maxIters);

relres = zeros(maxIters,1);
relres(1) = 1;

%% Solve for consistent ICs at k
if(Cur ~= Cur_prv)
    
    disp('Solving for consistent ICs');
    
    % Preallocate
    z_cons = zeros(Nz,maxIters);
    z_cons(:,1) = z;
    
    for idx = 1:(maxIters-1)
        
        % DAE eqns for current time-step
        [trash_var,g] = dae_dfn_zhouxinme2(x,z_cons(:,idx),Cur,p);   % g is the output values (concentrations/potentials/j_i/jsd)
        
        % Jacobian of DAE eqns
        [trash_var,trash_var,trash_var,g_z] = jac_dfn_zhouxinme2(x,z_cons(:,idx),Cur,p.f_x,p.f_z,p.g_x,p.g_z,p);
        
        % Newton Iteration
        Delta_z = -(g_z\g);   % Newton's method: x_n+1 = x_n - f(x_n)/f'(x_n); Delta_z = - f(x_n)/f'(x_n).
        z_cons(:,idx+1) = z_cons(:,idx) + Delta_z;
        
        % Check stopping criterion
        relres_z = norm(Delta_z,inf) / norm(z,inf);
        if(relres_z < tol)
            break;
        elseif(idx == (maxIters-1))
            fprintf(1,'Warning: Max Newton Iters Reached | RelChange = %3.2f%%\n',relres_z*100);
        end
        
    end
    
    z = z_cons(:,idx+1);
    
end

%% Solve Nonlinear System using Newton's Method
% DAE eqns for current time-step
[f,g] = dae_dfn_zhouxinme2(x,z,Cur,p);

% Initialize next x,z
x_nxt(:,1) = x;
z_nxt(:,1) = z;

% Iterate Newton's Method
FAC = [];
for idx = 1:(maxIters-1)

    % DAE eqns for next time-step
    [f_nxt, g_nxt] = dae_dfn_zhouxinme2(x_nxt(:,idx), z_nxt(:,idx), Cur_nxt, p);


    % Nonlinear system of equations
    F1 = x - x_nxt(:,idx) + p.delta_t/2 * (f + f_nxt);
    F2 = g_nxt;
    F = [F1; F2];

    % Jacobian of DAE eqns
    [f_x, f_z, g_x, g_z] = jac_dfn_zhouxinme2(x_nxt(:,idx),z_nxt(:,idx),Cur_nxt,p.f_x,p.f_z,p.g_x,p.g_z, p);
    
    % Jacobian of implicit time-stepped eqns
    F1_x = -speye(length(x)) + p.delta_t/2 * f_x;
    F1_z = p.delta_t/2 * f_z;
    F2_x = g_x;
    F2_z = g_z;
         
    J = [F1_x, F1_z; F2_x, F2_z];
      
    % Newman's method with y = [x; z]
    Delta_y = -(J\F); 

    x_nxt(:,idx+1) = x_nxt(:,idx) + Delta_y(1:Nx);
    z_nxt(:,idx+1) = z_nxt(:,idx) + Delta_y(Nx+1:end);
    
    % Check stopping criterion #1
    y = [x_nxt(:,idx+1); z_nxt(:,idx+1)];
    relres(idx+1) = norm(Delta_y,inf) / norm(y,inf);
    
    if( (relres(idx+1) < tol) && (norm(F,inf) < tol) )
        break;
    elseif(idx == (maxIters-1))
        fprintf(1,'Warning: Max Newton Iters Reached | RelChange = %1.4e%%\n',relres(end)*100);
    end
    
end

%% Output next x,z
x_nxtf = x_nxt(:,idx+1);
z_nxtf = z_nxt(:,idx+1);

newtonStats.iters = idx;
newtonStats.relres = relres;
newtonStats.condJac = condest(J);
varargout{1} = newtonStats;

