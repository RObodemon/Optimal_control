%% Main Function 
function midterm_dir_muti_shooting
close all;clear;clc;
tic % calculate computation time
% constant
u  = 1;
T  = 0.1405;
ve = 1.8758344;

nx      = 5; % # of states
K       = 16; % # of interval 2 4 8 16
n       = 2; % # degrees of polynomial 
A       = [];
B       = [];
Aeq     = [];
Beq     = [];
cguess  = zeros(K*(n+1),1);
t0      = 0;
tfguess = 3.5;

P_guess = ones(nx,K-1);
pmin    = -50*ones(nx,K-1);
pmax    = +50*ones(nx,K-1);
tfmin   = 0;
tfmax   = 20;
zmin    = [-50*ones(K*(n+1),1);pmin(:);tfmin];
zmax    = [+50*ones(K*(n+1),1);pmax(:);tfmax];

tau     = linspace(-1,+1,K+1);
zguess  = [cguess; P_guess(:);tfguess];
A       = [];
B       = [];
Aeq     = [];
Beq     = [];

option = optimset('Display','Iter','TolX',1e-8);
z      = fmincon(@Obj_Orbital,zguess,A,B,Aeq,Beq,zmin,zmax,@Error_Orbital,option,u,T,ve,K,n,nx,tau,t0)
toc
disp(['runtime: ',num2str(toc)]); % display computation time

[Eineq,Eeq,t,p] = Error_Orbital(z,u,T,ve,K,n,nx,tau,t0)
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
c    = z(1:end-1);
beta = polyval(c,t);
figure(3);
plot(t,beta,'ro-');
title('beta(t)');
legend('control angle');
grid on;
end

%% Error Function 
function [Eineq,Eeq,t,P] = Error_Orbital(z,u,T,ve,K,n,nx,tau,t0)

tf_G  = z(end);
c_G   = z(1:end-nx*(K-1)-1);
c_G   = reshape(c_G,n+1,K);
PTot0 = z(end-nx*(K-1):end-1);
PTot0 = reshape(PTot0,nx,K-1);
r0      = 1;               % initial position
rf      = 1.5;             % final positon
theta0  = 0;               % initial angle
vr0     = 0;               % initial linear velocity 
vrf     = 0;               % final linear velocity
vtheta0 = sqrt(u/r0);      % initial angular velocity  
vthetaf = sqrt(u/rf);      % final angular velocity
m0      = 1;               % initial mass
P0      = [r0; theta0; vr0; vtheta0; m0];

E = [];
t = [];
P = [];

options = odeset('RelTol',1e-8);
for k = 1:K
    if isequal(k,1)
        P0 = P0;
    else
        P0 = PTot0(:,k-1);
    end
    c           = c_G(:,k); % update polynomial
    tauspan     = [tau(k),tau(k+1)]; % update time span
    [tout,Pout] = ode113(@myOde,tauspan,P0,options,u,T,ve,c,K,nx,tau,t0,tf_G);
    ptf         = Pout(end,:).';
    if k<K
        E = [E;ptf-PTot0(:,k)];
    end
    t = [t; tout];
	P = [P; Pout];
end
r_tf      = ptf(1);
theta_tf  = ptf(2);
vr_tf     = ptf(3);
vtheta_tf = ptf(4);
m_tf      = ptf(5);
Eeq       = [E;r_tf-rf; vr_tf-vrf; vtheta_tf-vthetaf]; % collect all error
Eineq     = [];

end

%% Objectivd Function
function J = Obj_Orbital(z,u,T,ve,K,n,nx,tau,t0)

tf = z(end);
J  = tf;

end

%% Ode Function 
function Pdot = myOde(t,P,u,T,ve,c,K,nx,tau,t0,tf_G)

beta   = polyval(c,t); % control function
r      = P(1);
theta  = P(2);
vr     = P(3);
vtheta = P(4);
m      = P(5);

vtheta_dot = T*cos(beta)/m - vr*vtheta/r;
vr_dot     = vtheta^2/r + T*sin(beta)/m - u/(r^2);
theta_dot  = vtheta/r;
r_dot      = vr;
m_dot      = -T/ve;

Pdot = (tf_G-t0)/2*[r_dot; theta_dot; vr_dot; vtheta_dot; m_dot];

end