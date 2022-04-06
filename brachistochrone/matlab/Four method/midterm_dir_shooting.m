%% Main Function
function midterm_dir_shooting
close all;clear;clc;
tic % calculate computation time
% constant
u  = 1;
T  = 0.1405;
ve = 1.8758344;

n       = 5;         % # degrees of polynomial 
A       = [];
B       = [];
Aeq     = [];
Beq     = [];

cguess  = [0;0;0;0;1;0]; %zeros(n+1,1);% different guess of the coefficent will lead very different trajectory;
tfguess = 8;           %5,8 is the best
zguess  = [cguess; tfguess];
tfmin   = 0;
tfmax   = 20;
zmin    = [-100*ones(n+1,1);tfmin];
zmax    = [+100*ones(n+1,1);tfmax];

option  = optimset('Display','Iter','TolX',1e-8);
z       = fmincon(@Obj_Orbital,zguess,A,B,Aeq,Beq,zmin,zmax,@Error_Orbital,option,u,T,ve);
toc
disp(['runtime: ',num2str(toc)]); % display computation time

[Eineq,Eeq,t,p] = Error_Orbital(z,u,T,ve);
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
function [Eineq,Eeq,t,P] = Error_Orbital(z,u,T,ve)

tf = z(end);
c  = z(1:end-1);

r0      = 1;               % initial position
rf      = 1.5;             % final positon
theta0  = 0;               % initial angle
vr0     = 0;               % initial linear velocity 
vrf     = 0;               % final linear velocity
vtheta0 = sqrt(u/r0);      % initial angular velocity  
vthetaf = sqrt(u/rf);      % final angular velocity
m0      = 1;               % initial mass

P0 = [r0; theta0; vr0; vtheta0; m0];

tspan     = [0,tf];
options   = odeset('RelTol',1e-8);
[t,P]     = ode113(@myOde,tspan,P0,options,u,T,ve,c);
r_tf      = P(end,1);
theta_tf  = P(end,2);
vr_tf     = P(end,3);
vtheta_tf = P(end,4);
m_tf      = P(end,5);
Eeq       = [r_tf-rf; vr_tf-vrf; vtheta_tf-vthetaf;];
Eineq     = [];

end

%% Objectivd Function
function J = Obj_Orbital(z,u,T,ve)

tf = z(end);
J  = tf;

end

%% Ode Function 
function Pdot = myOde(t,P,u,T,ve,c)

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

Pdot = [r_dot; theta_dot; vr_dot; vtheta_dot; m_dot];

end