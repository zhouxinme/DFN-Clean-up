%% Reference Potential for Neg Electrode: Unref(theta)
%   Created July 12, 2011 by Scott Moura
% 
% Update Log by zhouxinme:
% Nov 5, 2014 - This reference potential function is taken from [Safari,
%               2011, JES], Modeling of a Commercial Graphite/LiFePO4 Cell.


function [Uref,varargout] = refPotentialAnode_LiFePO4_zhouxinme(p,theta)

if(~isreal(theta))
    beep;
    error('dfn:err','Complex theta');
%     pause;
end

% % Polynomail Fit
% Uref = ppvalFast(p.Uppn,theta);

% DUALFOIL: 0.01 < x < 0.9
a1 = 0.6379;
a2 = 0.5416;
a3 = -305.5309;
a4 = 0.044;
a5 = -0.1958;
a6 = -0.1088;
a7 = -0.1978;
a8 = -1.0571;
a9 = 0.0854;
a10 = -0.6875;
a11 = 0.0117;
a12 = 0.0529;
a13 = -0.0175;
a14 = -0.5692;
a15 = 0.0875;

Uref = a1 + a2.*exp(a3.*theta) ...
        + a4.*tanh((theta + a5)./a6) ...
        + a7.*tanh((theta + a8)./a9) ...
        + a10.*tanh((theta + a11)./a12) ...
        + a13.*tanh((theta + a14)./a15);

% Gradient of OCP wrt theta
if(nargout >= 2)

%     % Polynomial Fit
%     dUref = ppvalFast(p.dUppn,theta);
%     varargout{1} = dUref / p.c_s_n_max;


dUref = -165.47553544./exp(305.5309.*theta) - ...
        2.31615925058548.*sech(11.709601873536299.*(-1.0571 + theta)).^2 - ...
        0.2.*sech(11.428571428571429.*(-0.5692 + theta)).^2 - ...
        0.40441176470588236.*sech(9.191176470588236.*(-0.1958 + theta)).^2 - ...
        12.996219281663516.*sech(18.90359168241966.*(0.0117 + theta)).^2;

varargout{1} = dUref;

end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end

