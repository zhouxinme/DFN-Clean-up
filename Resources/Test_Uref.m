Ts = 0.01;

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

theta_p = 0:Ts:1;
theta_n = 1 - theta_p;

Uref_p = c1 + c2.*exp(c3.*(1 - theta_p).^c4) + c5.*exp(c6.*(1 - theta_p).^c7) ...
        + c8.*exp(c9.*(1 - theta_p).^c10);
    
dUref_p = (-89.26349843079201.*(1 - theta_p).^0.3198000000000001)./exp(80.2493.*(1 - ...
        theta_p).^1.3198) - 0.00025009628839914.*exp(20.2646.*(1 - ...
        theta_p).^3.7995).*(1 - theta_p).^2.7995 + ...
        0.00025008610382119.*exp(20.2645.*(1 - theta_p).^3.8003).*(1 - ...
        theta_p).^2.8003;

Uref_p_1 = zeros(size(Uref_p));
Uref_p_1(:,1) = Uref_p(1);
Uref_p_1(:,2:end) = Uref_p(:,1:end-1)+ Ts*dUref_p(:,1:end-1);

Uref_p_2 = zeros(size(Uref_p));
Uref_p_2(:,1) = Uref_p(1);
Uref_p_2(:,2:end) = Uref_p(:,1:end-1) + 1/2*Ts*(dUref_p(:,1:end-1)+dUref_p(:,2:end));

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

Uref_n = a1 + a2.*exp(a3.*theta_n) ...
        + a4.*tanh((theta_n + a5)./a6) ...
        + a7.*tanh((theta_n + a8)./a9) ...
        + a10.*tanh((theta_n + a11)./a12) ...
        + a13.*tanh((theta_n + a14)./a15);
    
dUref_n = -165.47553544./exp(305.5309.*theta_n) - ...
        2.31615925058548.*sech(11.709601873536299.*(-1.0571 + theta_n)).^2 - ...
        0.2.*sech(11.428571428571429.*(-0.5692 + theta_n)).^2 - ...
        0.40441176470588236.*sech(9.191176470588236.*(-0.1958 + theta_n)).^2 - ...
        12.996219281663516.*sech(18.90359168241966.*(0.0117 + theta_n)).^2;
    
Uref_n_1 = zeros(size(Uref_n));
Uref_n_1(:,1) = Uref_n(1);
Uref_n_1(:,2:end) = Uref_n(:,1:end-1)+ Ts*dUref_n(:,1:end-1);

Uref_n_2 = zeros(size(Uref_n));
Uref_n_2(:,1) = Uref_n(1);
Uref_n_2(:,2:end) = Uref_n(:,1:end-1) + 1/2*Ts*(dUref_n(:,1:end-1)+dUref_n(:,2:end));

delta_Uref = Uref_p - Uref_n;

figure
plot(theta_p,Uref_p,'b', ...
    theta_p,Uref_p_1,'k', ...
     theta_p,Uref_p_2,'r', ...
    'linewidth',1.5)
hold on
plot(theta_p,Uref_n,'b--', ...
    theta_p,Uref_n_1,'k--', ...
    theta_p,Uref_n_2,'r--', ...
    'linewidth',1.5)

figure
plot(theta_n,delta_Uref,'linewidth',1.5)

figure
plot(theta_p,dUref_p,'b', ...
    theta_p,dUref_n,'r', ...
    'linewidth',1.5)