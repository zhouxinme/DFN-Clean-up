%% Reference Potential for Pos. Electrode: Unref(theta_n)
%   Created July 12, 2011 by Scott Moura
% 
% Update Log by zhouxinme:
% Nov 5, 2014 - This reference potential function is taken from [Safari,
%               2011, JES], Modeling of a Commercial Graphite/LiFePO4 Cell.

function [Uref,varargout] = refPotentialCathode_LiFePO4_zhouxinme(p,theta)

% % Polynomail Fit
% Uref = ppvalFast(p.Uppp,theta);

% DUALFOIL: LFP
c1 = 3.4323;
c2 = -0.8428;
c3 = -80.2493;
c4 = 1.3198;
c5 = -3.2474e-6;
c6 = 20.2645;
c7 = 3.8003;
c8 = 3.2482e-6;
c9 = 20.2646;
c10 = 3.7995;

Uref = c1 + c2.*exp(c3.*(1 - theta).^c4) + c5.*exp(c6.*(1 - theta).^c7) ...
        + c8.*exp(c9.*(1 - theta).^c10);

% Gradient of OCP wrt theta
if(nargout == 2)
    
%     % Polynomail Fit
%     dUref = ppvalFast(p.dUppp,theta);
%     varargout{1} = dUref / p.c_s_p_max;

dUref = (-89.26349843079201.*(1 - theta).^0.3198000000000001)./exp(80.2493.*(1 - ...
        theta).^1.3198) - 0.00025009628839914.*exp(20.2646.*(1 - ...
        theta).^3.7995).*(1 - theta).^2.7995 + ...
        0.00025008610382119.*exp(20.2645.*(1 - theta).^3.8003).*(1 - ...
        theta).^2.8003;

varargout{1} = dUref;
    
end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end