%% Main Function
function midterm_ind_muti_shooting
close all;clear;clc
tic % calculate computation time
% constant
u  = 1;
T  = 0.1405;
ve = 1.8758344;

r0      = 1;               % initial position
rf      = 1.5;             % final positon
theta0  = 0;               % initial angle
vr0     = 0;               % initial linear velocity 
vrf     = 0;               % final linear velocity
vtheta0 = sqrt(u/r0);      % initial angular velocity  
vthetaf = sqrt(u/rf);      % final angular velocity
m0      = 1;               % initial mass
nx      = 5;               % # of states;

K    = 16;  % 2,4,8,16
tau  = linspace(-1,+1,K+1);
t0   = 0;
tf_G = 3.5; % the guess for the time is very important,can not bigger than 11
lamG = [-1 0 0 2 0]; % the first and the forth guess is very important   

lam_r0_G       = lamG(1);
lam_theta0_G   = lamG(2);
lam_vr0_G      = lamG(3);
lam_vtheta0_G  = lamG(4);
lam_m_G        = lamG(5);

lam_G = [lam_r0_G; lam_theta0_G; lam_vr0_G; lam_vtheta0_G; lam_m_G; tf_G];
P0_G  = ones(2*nx,K-1);
z_G   = [lam_G;P0_G(:)];

options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
z       = fsolve(@Error_Orbital,z_G,options,r0,rf,theta0,vr0,vrf,vtheta0,vthetaf,m0,t0,u,T,ve,nx,tau,K)

toc
disp(['runtime: ',num2str(toc)]); % display computation time

[E,t,p] = Error_Orbital(z,r0,rf,theta0,vr0,vrf,vtheta0,vthetaf,m0,t0,u,T,ve,nx,tau,K)
% plot
figure(1);
plot(t,p(:,1),'r*-');
hold on;
plot(t,p(:,2),'go-');
hold on;
plot(t,p(:,3),'b*-');
hold on;
plot(t,p(:,4),'mo-');
hold on;
plot(t,p(:,5),'k+-');
grid on;
legend('r','theta','vr','vtheta','m');
title('states(t)');

% plot trajectory
figure(2);
polarplot(p(:,2),p(:,1),'r');
title('trajectory');
hold on;

% plot control
figure(3);
beta    = atan2(p(:,8),p(:,9));
beta    = unwrap(beta);
plot(t,beta,'ro-');
title('beta(t)');
grid on;
end

%% Error Function
function [E,t,P] = Error_Orbital(z,r0,rf,theta0,vr0,vrf,vtheta0,vthetaf,m0,t0,u,T,ve,nx,tau,K)

lam_r0_G      = z(1);
lam_theta0_G  = z(2);
lam_vr0_G     = z(3); 
lam_vtheta0_G = z(4); 
lam_m0_G      = z(5); 
tf_G          = z(6);

P0    = [r0; theta0; vr0; vtheta0; m0; lam_r0_G; lam_theta0_G; lam_vr0_G; lam_vtheta0_G; lam_m0_G];
PTot0 = z(7:end);
PTot0 = reshape(PTot0,2*nx,K-1);

E = [];
t = [];
P = [];
option = odeset('RelTol',1e-8);

for k = 1:K
    if isequal(k,1);
        P0 = P0;
    else
        P0 = PTot0(:,k-1);
    end

    tauspan      = [tau(k), tau(k+1)];
    [tout, Pout] = ode113(@myOde_Orbital,tauspan,P0,option,u,T,ve,t0,tf_G);
    Ptf          = Pout(end,:).';

    if k<K
        E = [E;Ptf-PTot0(:,k)];
    end
    t = [t;tout];
    P = [P;Pout];
end

r_tf      = Ptf(1);
theta_tf  = Ptf(2);
vr_tf     = Ptf(3);
vtheta_tf = Ptf(4);
m_tf      = Ptf(5);

lam_r_tf      = Ptf(6);
lam_theta_tf  = Ptf(7);
lam_vr_tf     = Ptf(8);
lam_vtheta_tf = Ptf(9);
lam_m_tf      = Ptf(10);

beta    = atan2(lam_vr_tf,lam_vtheta_tf);
beta    = unwrap(beta);

Htf = lam_r_tf*vr_tf +...
      lam_theta_tf*vtheta_tf/r_tf +...
      lam_vr_tf*(vtheta_tf^2/r_tf + T*sin(beta)/m_tf - u/r_tf^2) +...
      lam_vtheta_tf*(T*cos(beta)/m_tf - vr_tf*vtheta_tf/r_tf) +...
      lam_m_tf*(-T/ve);

E   = [E;r_tf-rf; vr_tf-vrf; vtheta_tf-vthetaf; lam_theta_tf; lam_m_tf-1; Htf];

end

%% state and costate Function
function P_dot = myOde_Orbital(tspan,P,u,T,ve,t0,tf)

r          = P(1);
vr         = P(3);
vtheta     = P(4);
m          = P(5);
lam_r      = P(6);
lam_theta  = P(7);
lam_vr     = P(8);
lam_vtheta = P(9);

beta    = atan2(lam_vr,lam_vtheta);
beta    = unwrap(beta);

r_dot      = vr;
theta_dot  = vtheta/r;
vr_dot     = vtheta^2/r + T*sin(beta)/m - u/(r^2);
vtheta_dot = T*cos(beta)/m - vr*vtheta/r;
m_dot      = -T/ve;

lam_r_dot      = -lam_vtheta*vr*vtheta/r^2 + lam_vr*(vtheta^2/r^2 - 2*u/r^3) + lam_theta*vtheta/r^2;
lam_theta_dot  = 0;
lam_vr_dot     = lam_vtheta*vtheta/r - lam_r;
lam_vtheta_dot = lam_vtheta*vr/r - 2*lam_vr*vtheta/r - lam_theta/r;
lam_m_dot      = lam_vtheta*T*cos(beta)/m^2 + lam_vr*T*sin(beta)/m^2;

P_dot = (tf-t0)/2*[r_dot; theta_dot; vr_dot; vtheta_dot; m_dot;...
         lam_r_dot; lam_theta_dot; lam_vr_dot; lam_vtheta_dot; lam_m_dot];
end
