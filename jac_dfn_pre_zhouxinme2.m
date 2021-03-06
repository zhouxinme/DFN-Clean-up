%% Jacobian for Doyle-Fuller-Newman Model
%   Created May 30, 2012 by Scott Moura
%   Preprocessed elements, i.e. non-state dependent

function [f_x, f_z, g_x, g_z] = jac_dfn_pre_zhouxinme2(p)


%% Parse out states

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

ind_csn = 1:Ncsn;
ind_csp = Ncsn+1:Ncsn+Ncsp;

ind_ce = Ncsn+Ncsp+1:Nc;

ind_phi_s_n = 1:Nn;
ind_phi_s_p = Nn+1:Nnp;

ind_ien = Nnp+1:Nnp+Nn;
ind_iep = Nnp+Nn+1:2*Nnp;

ind_phi_e = 2*Nnp+1 : 2*Nnp+Nce;

ind_jp = 2*Nnp+Nce+1 : 2*Nnp+Nce+Np;

ind_jn1 = 2*Nnp+Nce+Np+1 : 3*Nnp+Nce;
ind_jsd = 3*Nnp+Nce+1 : 3*Nnp+Nce+Nn;

%% Preallocate Jacobian
f_x = zeros(Nx);
f_z = zeros(Nx,Nz);
g_x = zeros(Nz,Nx);
g_z = zeros(Nz);

%% Li Diffusion in Solid Phase: c_s(x,r,t)
% Preallocate
Cell_Acsn = cell(Nn,1);
Cell_Bcsn = cell(Nn,1);
Cell_Acsp = cell(Np,1);
Cell_Bcsp = cell(Np,1);

% Loop through each "comb tooth" in anode
for idx = 1:Nn
    Cell_Acsn{idx} = p.A_csn;
    Cell_Bcsn{idx} = p.B_csn;
end
f_x(ind_csn,ind_csn) = blkdiag(Cell_Acsn{:});
f_z(ind_csn,ind_jn1) = blkdiag(Cell_Bcsn{:});   % csn is only affected by jn1, no effects of jsd.

% Loop through each "comb tooth" in cathode
for idx = 1:Np
    Cell_Acsp{idx} = p.A_csp;
    Cell_Bcsp{idx} = p.B_csp;
end
f_x(ind_csp,ind_csp) = blkdiag(Cell_Acsp{:});
f_z(ind_csp,ind_jp) = blkdiag(Cell_Bcsp{:});

%% Potential in Solid Phase: phi_s(x,t)
g_z(ind_phi_s_n,ind_phi_s_n) = p.F1_psn;
g_z(ind_phi_s_p,ind_phi_s_p) = p.F1_psp;

g_z(ind_phi_s_n,ind_ien) = p.F2_psn(:,2:end-1);
g_z(ind_phi_s_p,ind_iep) = p.F2_psp(:,2:end-1);

%% Electrolyte Current: i_e(x,t)
g_z(ind_ien,ind_ien) = p.F1_ien;
g_z(ind_iep,ind_iep) = p.F1_iep;

g_z(ind_ien,ind_jn1) = p.F2_ien;
g_z(ind_ien,ind_jsd) = p.F2_ien;
g_z(ind_iep,ind_jp) = p.F2_iep;

%%
f_x = sparse(f_x);
f_z = sparse(f_z);
g_x = sparse(g_x);
g_z = sparse(g_z);

% spy(Jac)
% pause;
