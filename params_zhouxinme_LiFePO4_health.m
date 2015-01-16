%% Parameters for Electrochemical Model from the Literature
%   Xin Zhou update from SetParav20.m and para in Joel's code.
% 
% Update Log by zhouxinme:
% Nov 4th 2014 - Change all the parameters to what is inside SetPara. The
%                name in SetParav20 and in Joel's code is also included in
%                the comments.

%% Geometric Params
% Thickness of each layer
p.L_n = 2.8853e-5;     % Thickness of negative electrode [m]   % delta_n   %lengthsec
p.L_s = 1.6971e-5;     % Thickness of separator [m]   % delta_p
p.L_p = 6.5205e-5;     % Thickness of positive electrode [m]   % delta_sep

L_ccn = 25e-6;    % Thickness of negative current collector [m]   % NOT in SetPara   % NOT in para
L_ccp = 25e-6;    % Thickness of negative current collector [m]   % NOT in SetPara

% Particle Radii
p.R_s_n = 3.5961e-6;   % Radius of solid particles in negative electrode [m]   % R_n   % SphericalRadius
p.R_s_p = 1.6371e-7;   % Radius of solid particles in positive electrode [m]   % R_p

% Volume fractions
p.epsilon_s_n = 0.381;      % Volume fraction in solid for neg. electrode   % epsilon_s_n   % epsilon1
p.epsilon_s_p = 0.48;      % Volume fraction in solid for pos. electrode   % epsilon_s_p

p.epsilon_e_n = 0.619;   % Volume fraction in electrolyte for neg. electrode   % epsilon_e_n   % epsilon2
p.epsilon_e_s = 0.3041;   % Volume fraction in electrolyte for separator   % epsilon_e_sep
p.epsilon_e_p = 0.52;   % Volume fraction in electrolyte for pos. electrode   % epsilon_e_p

epsilon_f_n = 0.1;  % Volume fraction of filler in neg. electrode % NOT in SetPara   % NOT in para
epsilon_f_p = 0.2;  % Volume fraction of filler in pos. electrode % NOT in SetPara

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]   % a_s_n
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]   % a_s_p

% Mass densities
rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]   % NOT in SetPara   % NOT in para
rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]   % NOT in SetPara   % NOT in para
rho_e =  1324;    % Electrolyte [kg/m^3]   % NOT in SetPara   % NOT in para
rho_f = 1800;     % Filler [kg/m^3]   % NOT in SetPara   % NOT in para
rho_ccn = 8954;   % Current collector in negative electrode   % NOT in SetPara   % NOT in para
rho_ccp = 2707;   % Current collector in positive electrode   % NOT in SetPara   % NOT in para

% Compute cell mass [kg/m^2]
m_n = p.L_n * (rho_e*p.epsilon_e_n + rho_sn*p.epsilon_s_n + rho_f*epsilon_f_n);
m_s = p.L_s * (rho_e*p.epsilon_e_n);
m_p = p.L_p * (rho_e*p.epsilon_e_p + rho_sp*p.epsilon_s_p + rho_f*epsilon_f_p);
m_cc = rho_ccn*L_ccn + rho_ccp*L_ccp;

% Lumped density [kg/m^2]
p.rho_avg = m_n + m_s + m_p + m_cc;
% p.rho_avg = 1e6; % (rho_avg * c_P lumped into rho_avg)   % NOT in SetPara   % NOT in para

%% Transport Params
% Diffusion coefficient in solid
p.D_s_n = 8.2557e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]   % D_s_n   % d1
p.D_s_p = 1.7362e-14;  % Diffusion coeff for solid in pos. electrode, [m^2/s]   % D_s_p

% Diffusion coefficient in electrolyte
p.D_e = 6.9114e-10;    % Diffusion coeff for electrolyte, [m^2/s]   % NOT in SetPara   % d2

p.brug = 1.452;       % Bruggeman porosity   % NOT in SetPara   % brug

% Conductivity of solid
p.sig_n = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]   % NOT in SetPara   % sigma
p.sig_p = 100;    % Conductivity of solid in pos. electrode, [1/Ohms*m]

p.sig_eff_n = p.sig_n * p.epsilon_s_n^p.brug;    % Eff. conductivity in neg. electrode, [1/Ohms*m]
p.sig_eff_p = p.sig_p * p.epsilon_s_p^p.brug;    % Eff. conductivity in pos. electrode, [1/Ohms*m]

% Conductivity of electrolyte

% Miscellaneous
p.t_plus = 0.2495;       % Transference number   % Not in SetPara   % NOT in SetPara   % transference
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]   % F   % Faraday
p.Area = 0.3108;           % Electrode current collector area [m^2]   % A   % Joel's code - SISO_subModelId: vol = 0.013^2*pi*0.065; A = vol/(out1(1) + out1(2) + out1(3))

%% Kinetic Params
p.R = 8.314;       % Gas constant, [J/mol-K]   % R   % UniversalGasConstant
p.alph = 0.5;         % Charge transfer coefficients   % alpha
p.R_f_n = 0.0037;       % Resistivity of SEI layer, [Ohms*m^2]   % R_SEI   % rSEI
p.R_f_p = 0;       % Resistivity of SEI layer, [Ohms*m^2]

% Reaction rates
p.k_n = 8.6963e-7;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]   % k_n   % krate
p.k_p = 1.1267e-7; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]   % k_p

%% Thermodynamic Params
% Cubic Splines of Equilibrium potentials (ws/pp_LiFEPO4.mat)
% load('ws/BoschOCP');
% 
% ppn = spline(theta,OCPn);
% dppn = ppn;
% dppn.coefs = [3*ppn.coefs(:,1), 2*ppn.coefs(:,2), ppn.coefs(:,3)];
% dppn.order = 3;
% 
% ppp = spline(theta,OCPp);
% dppp = ppp;
% dppp.coefs = [3*ppp.coefs(:,1), 2*ppp.coefs(:,2), ppp.coefs(:,3)];
% dppp.order = 3;
% 
% p.Uppn = ppn;
% p.dUppn = dppn;
% p.Uppp = ppp;
% p.dUppp = dppp;

% Thermal dynamics
p.C_p = 1;   % Heat capacity, [J/kg-K] (rho_avg * c_P lumped into rho_avg)   % NOT in SetPara   % NOT in para
p.h = 20;   % Heat transfer coefficient, [W/K-m^2]   % NOT in SetPara   % NOT in para

% Ambient Temperature
p.T_amp = 298.15; % [K]   % T   % Temperature

% Entropy coefficients
p.dUref_dT = -0.4e-3; % [V/K] approx. from Al Hallaj et al 2000, JPS   % NOT in SetPara   % NOT in para

%% Concentrations
p.c_s_n_max = 2.9482e4;   % Max concentration in anode, [mol/m^3]   % c_s_n_max   % c1max
p.c_s_p_max = 1.0355e4;    % Max concentration in cathode, [mol/m^3]   % c_s_p_max

p.n_Li_s = [1/2*(p.epsilon_s_n*p.L_n*p.c_s_n_max + p.epsilon_s_p*p.L_p*p.c_s_p_max)]*p.Area;        % Total moles of lithium in solid phase [mol]   % I just guess this value since it is not available
p.c_e = 1.2669e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]   % c_e_true

%% Cutoff voltages
p.volt_max = 3.6;
p.volt_min = 2.0;

%% Discretization parameters
% Discrete time step
p.delta_t = 1e-2;

% Pade Order
p.PadeOrder = 3;

% Finite difference points along r-coordinate
% p.Nr = 10;
% p.delta_r_n = p.R_s_n / p.Nr;
% p.delta_r_p = p.R_s_p / p.Nr;

% Finite difference points along x-coordinate
p.Nxn = 10;
p.Nxs = 5;
p.Nxp = 10;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

%% Health
p.i_0sd = 1.5e-6;   % Exchange current density, [A/m^2], Ramadass: 1.5e-6; Joel: 4e-8
p.Urefsd = 0.4;   % Reference potential of the side reaction, [V], from Ramadass et al 2004, JES
