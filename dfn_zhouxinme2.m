%% Doyle-Fuller-Newman Model
%   Created May 22, 2012 by Scott Moura; Edited Jan 15th, 2014 by Xin Zhou.
% 
% State vectors:
%  x0 = [c_s_n(:,1); c_s_p(:,1); c_e(:,1); T(1)];
%  z0 = [phi_s_n(:,1); phi_s_p(:,1); i_en(:,1); i_ep(:,1);...
%        phi_e(:,1); jp(:,1); jn1(:,1); jsd(:,1)];
% 
% Update Log by zhouxinme:
% Jan 15th 2014: Clean up the DFN model.

clc;
clear;
tic;

%% Model Construction
% Electrochemical Model Parameters
run params_zhouxinme_LiFePO4_health

% Vector lengths
Ncsn = p.PadeOrder * (p.Nxn-1);   % Length of c_{s,n}, (n_r)*(n_x).
Ncsp = p.PadeOrder * (p.Nxp-1);   % Length of c_{s,p}, (n_r)*(n_x).
Nce = p.Nx - 3;   % Length of c_e.
Nc = Ncsn+Ncsp+Nce;   % Length of all concentration states.
Nn = p.Nxn - 1;   % # of nodes in the anode.
Np = p.Nxp - 1;   % # of nodes in the cathode.
Nnp = Nn+Np;   % # of nodes in the electrodes.

Nx = Ncsn + Ncsp + Nce + 1;   % # of states in x.
Nz = 3*Nnp + Nce + Nn;   % # of states in z.

%% Input Signal

% Choose the input signal type
IST = 1;

switch IST
    case 1   % Manual Data
        t = -0.05:p.delta_t:1;       % 60*60
        Iamp = zeros(length(t),1);
        Iamp(t >= 0) = -20*2.3;
        % Iamp(t >= (20*60)) = 0;
        % Iamp(t >= 20) = 0;
        % Iamp(t >= 30) = 10;
        % Iamp(t >= 40) = 5;
        I = Iamp;
    case 2   % Pulse Data
        t = -2:p.delta_t:240;
        I(t >= 0) = 0;
        I(mod(t,40) < 20) = 175;
        Iamp = I;
    case 3   % Experimental Data
         load('data/UDDSx2_batt_ObsData.mat');
         tdata = t;
         Tfinal = tdata(end);
         t = -2:p.delta_t:Tfinal;
         Iamp = interp1(tdata,I,t,'spline',0);
         Ah_amp = trapz(tdata,I)/3600;
         I = Iamp * (0.4*35)/Ah_amp;
         %cut the simulation time to 500 seconds
         t=-2:p.delta_t:500;
         I=I(1:length(t));
end

NT = length(t);

%% Initial Conditions & Preallocation
% Solid concentration
V0 = 3.5;
[csn0,csp0] = init_cs_zhouxinme(p,V0);

c_s_n0 = zeros(p.PadeOrder,1);
c_s_p0 = zeros(p.PadeOrder,1);

%%%%% Initial condition based on controllable canonical form
% c_s_n0(1) = csn0 * (-p.R_s_n/3) * (p.R_s_n^4 / (3465 * p.D_s_n^2));
% c_s_p0(1) = csp0 * (-p.R_s_p/3) * (p.R_s_p^4 / (3465 * p.D_s_p^2));

%%%%% Initial condition based on Jordan form
c_s_n0(3) = csn0;   % c_s_n0(3): for 3rd order Pade approximation
c_s_p0(3) = csp0;   % c_s_p0(3): for 3rd order Pade approximation
%%%%%

c_s_n = zeros(Ncsn,NT);
c_s_p = zeros(Ncsp,NT);

c_s_n(:,1) = repmat(c_s_n0, [Nn 1]);
c_s_p(:,1) = repmat(c_s_p0, [Np 1]);

% Electrolyte concentration
c_e = zeros(Nce,NT);
c_e(:,1) = p.c_e * ones(Nce,1);

c_ex = zeros(Nce+4,NT);   % # of all the nodes in x direction is (p.Nx+1)
c_ex(:,1) = c_e(1,1) * ones(Nce+4,1);

% Temperature
T = zeros(1,NT);   % Uniform distribution of temperature
T(1) = p.T_amp;

% Solid Potential
Uref_n0 = refPotentialAnode_LiFePO4_zhouxinme(p, csn0(1)*ones(Nn,1) / p.c_s_n_max);
Uref_p0 = refPotentialCathode_LiFePO4_zhouxinme(p, csp0(1)*ones(Np,1) / p.c_s_p_max);

phi_s_n = zeros(Nn,NT);
phi_s_p = zeros(Np,NT);
phi_s_n(:,1) = Uref_n0;
phi_s_p(:,1) = Uref_p0;

% Electrolyte Current
i_en = zeros(Nn,NT);
i_ep = zeros(Np,NT);

% Electrolyte Potential
phi_e = zeros(Nce,NT);

% Molar Ionic Flux in the cathode
jp = zeros(Np,NT);

% Health: jn1 + jsd = jn; hence do not need an additional state for molar ionic flux in the anode
jn1 = zeros(Nn,NT);
jsd = zeros(Nn,NT);

% Surface concentration
c_ss_n = zeros(Nn,NT);
c_ss_p = zeros(Np,NT);
c_ss_n(:,1) = repmat(csn0, [Nn 1]);
c_ss_p(:,1) = repmat(csp0, [Np 1]);

% Volume average concentration
c_avg_n = zeros(Nn,NT);
c_avg_p = zeros(Np,NT);
c_avg_n(:,1) = repmat(csn0, [Nn 1]);
c_avg_p(:,1) = repmat(csp0, [Np 1]);

SOC = zeros(1,NT);
SOC(1) = mean(c_avg_n(:,1)) / p.c_s_n_max;

% Overpotential
eta_n = zeros(Nn,NT);
eta_p = zeros(Np,NT);

% Constraint Outputs
c_e_0p = zeros(1,NT);
c_e_0p(1) = c_ex(1,1);

eta_s_Ln = zeros(1,NT);
eta_s_Ln(1) = phi_s_p(1,1) - phi_e(1,1);

% Voltage
Volt = zeros(1,NT);
Volt(1) = phi_s_p(end,1) - phi_s_n(1,1);

% Conservation of Li-ion matter
nLi = zeros(1,NT);
nLidot = zeros(1,NT);

% Stats (Newman method)
newtonStats.iters = zeros(NT,1);
newtonStats.relres = cell(NT,1);
newtonStats.condJac = zeros(NT,1);

% Initial Conditions
x0 = [c_s_n(:,1); c_s_p(:,1); c_e(:,1); T(1)];

z0 = [phi_s_n(:,1); phi_s_p(:,1); i_en(:,1); i_ep(:,1);...
      phi_e(:,1); jp(:,1); jn1(:,1); jsd(:,1)];

%% Preallocate
x = zeros(Nx, NT);
z = zeros(Nz, NT);

x(:,1) = x0;
z(:,1) = z0;

%% Precompute data
% Solid concentration matrices
[A_csn,B_csn,A_csp,B_csp,C_csn,C_csp,A_csn_normalized, A_csp_normalized] = c_s_mats(p);
p.A_csn = A_csn;
p.A_csn_normalized= A_csn_normalized;
p.B_csn = B_csn;
p.A_csp = A_csp;
p.A_csp_normalized=A_csp_normalized;
p.B_csp = B_csp;
p.C_csn = C_csn;
p.C_csp = C_csp;

% Electrolyte concentration matrices
[trash_var,trash_var,C_ce] = c_e_mats_federico(p,c_ex);
p.C_ce = C_ce;

% Solid Potential
[F1_psn,F1_psp,F2_psn,F2_psp,G_psn,G_psp,...
    C_psn,C_psp,D_psn,D_psp] = phi_s_mats(p);
p.F1_psn = F1_psn;
p.F1_psp = F1_psp;
p.F2_psn = F2_psn;
p.F2_psp = F2_psp;
p.G_psn = G_psn;
p.G_psp = G_psp;
p.C_psn = C_psn;
p.C_psp = C_psp;
p.D_psn = D_psn;
p.D_psp = D_psp;

% Electrolyte Current
[F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p);
p.F1_ien = F1_ien;
p.F1_iep = F1_iep;
p.F2_ien = F2_ien;
p.F2_iep = F2_iep;
p.F3_ien = F3_ien;
p.F3_iep = F3_iep;

% Jacobian
[f_x, f_z, g_x, g_z] = jac_dfn_pre_zhouxinme2(p);
p.f_x = f_x;
p.f_z = f_z;
p.g_x = g_x;
p.g_z = g_z;
clear f_x f_z g_x g_z

% Calculate more states
jn = zeros(Nn,NT);

% Average current densities
Jsd = zeros(1,NT);
Jn1 = zeros(1,NT);
Jn = zeros(1,NT);
Jp = zeros(1,NT);

%% Integrate!
disp('Simulating DFN Model...');

for k = 1:(NT-1)
    
    % Current
    if(k == 1)
        Cur_vec = [I(k), I(k), I(k+1)];
    else
        Cur_vec = [I(k-1), I(k), I(k+1)];
    end
    
    % Step-forward in time
    [x(:,k+1), z(:,k+1), stats] = cn_dfn_zhouxinme2(x(:,k),z(:,k),Cur_vec,p);

    % Parse out States
    c_s_n(:,k+1) = x(1:Ncsn, k+1);
    c_s_p(:,k+1) = x(Ncsn+1:Ncsn+Ncsp, k+1);
    c_e(:,k+1) = x(Ncsn+Ncsp+1:Nc, k+1);
    T(k+1) = x(end, k+1);
    
    phi_s_n(:,k+1) = z(1:Nn, k+1);
    phi_s_p(:,k+1) = z(Nn+1:Nnp, k+1);
    i_en(:,k+1) = z(Nnp+1:Nnp+Nn, k+1);
    i_ep(:,k+1) = z(Nnp+Nn+1:2*Nnp, k+1);
    phi_e(:,k+1) = z(2*Nnp+1:2*Nnp+Nce, k+1);

    jp(:,k+1) = z(2*Nnp+Nce+1:2*Nnp+Nce+Np, k+1);
    jn1(:,k+1) = z(2*Nnp+Nce+Np+1:3*Nnp+Nce, k+1);
    jsd(:,k+1) = z(3*Nnp+Nce+1:3*Nnp+Nce+Nn, k+1);
    
    jn(:,k+1) = jn1(:,k+1) + jsd(:,k+1);
    
    Jsd(k+1) = sum(jsd(:,k+1))*p.Faraday*p.a_s_n/Nn;
    Jn1(k+1) = sum(jn1(:,k+1))*p.Faraday*p.a_s_n/Nn;
    Jn(k+1) = sum(jn(:,k+1))*p.Faraday*p.a_s_n/Nn;
    Jp(k+1) = sum(jp(:,k+1))*p.Faraday*p.a_s_p/Np;
    
    newtonStats.iters(k+1) = stats.iters;
    newtonStats.relres{k+1} = stats.relres;
    newtonStats.condJac(k+1) = stats.condJac;
    
    % Output data
    [trash_var, trash_var, y] = dae_dfn_zhouxinme2(x(:,k+1),z(:,k+1),I(k+1),p);
    
    c_ss_n(:,k+1) = y(1:Nn);
    c_ss_p(:,k+1) = y(Nn+1:Nnp);
    
    c_avg_n(:,k+1) = y(Nnp+1:Nnp+Nn);
    c_avg_p(:,k+1) = y(Nnp+Nn+1 : 2*Nnp);
    SOC(k+1) = mean(c_avg_n(:,k+1)) / p.c_s_n_max;
    
    c_ex(:,k+1) = y(2*Nnp+1:2*Nnp+Nce+4);
    
    eta_n(:,k+1) = y(2*Nnp+Nce+4+1 : 2*Nnp+Nce+4+Nn);
    eta_p(:,k+1) = y(2*Nnp+Nce+4+Nn+1 : 3*Nnp+Nce+4);
    
    c_e_0p(k+1) = y(3*Nnp+Nce+4+2);
    eta_s_Ln(k+1) = y(3*Nnp+Nce+4+3);
    
    Volt(k+1) = y(3*Nnp+Nce+4+4);
    nLi(k+1) = y(3*Nnp+Nce+4+5);
    nLidot(k+1) = y(3*Nnp+Nce+4+6);
    
    eta_s_n = phi_s_n - phi_e(1:Nn,:);
    eta_s_p = phi_s_p - phi_e(end-Np+1:end, :);
    
    fprintf(1,'Time : %3.2f sec | Current : %2.4f A/m^2 | SOC : %1.3f | Voltage : %2.4fV\n',...
        t(k),I(k+1),SOC(k+1),Volt(k+1));
    
    if(Volt(k+1) < p.volt_min)
        fprintf(1,'Min Voltage of %1.1fV exceeded\n',p.volt_min);
        beep;
        break;
    elseif(Volt(k+1) > p.volt_max)
        fprintf(1,'Max Voltage of %1.1fV exceeded\n',p.volt_max);
        beep;
        break;
    elseif(any(c_ex(:,k) < 1))
        fprintf(1,'c_e depleted below 1 mol/m^3\n');
        beep;
        break;
    end

end


%% Outputs
disp('Simulating Output Vars...');
simTime = toc;
fprintf(1,'Simulation Time : %3.2f min\n',simTime/60);
