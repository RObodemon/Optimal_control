%% main function
function Brachistochrone_OCP_main

nx = 3; % # of states
x0 = 0;
xf = 2;
y0 = 0;
yf = 2;
v0 = 0;
t0 = 0;
lamx_G = 0;
lamy_G = 0;
lamv_G = 0;
tf_G = 10;
g = 10;
K = 2;
tau = linspace(-1,+1,K+1);
P0_G = zeros(2*nx,K-1);
lam_G = [lamx_G, lamy_G, lamv_G,tf_G]';
z_G = [lam_G;P0_G(:)];

options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
z = fsolve(@Error_Brac,z_G,options,t0,x0,y0,v0,xf,yf,g,nx,K,tau)
[E,t,p] = Error_Brac(z,t0,x0,y0,v0,xf,yf,g,nx,K,tau)


% plot all 
figure(1);
plot(t,p(:,1),'r*-');   % x(t)
hold on;
plot(t,p(:,2),'bo-');   % y(t)
hold on;
plot(t,p(:,3),'go-');   % v(t)
hold on;
grid on;
legend('x(t)','y(t)','v(t)');

% plot trajectory
figure(2);
plot(p(:,1),p(:,2),'ro-');
title('Brachistochrone y/x');
grid;

% plot theta(t)
theta_G = 0;
theta = [];
options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
for i = 1:length(t)
    theta(i) = fsolve(@solveControl,theta_G,options,p(i,4),p(i,5),p(i,6),p(i,3),g);
end
figure(3);
plot(t,theta,'bo-')
title('theta(t)');
grid on;

end



%% Error function
function [E,t,p] = Error_Brac(z,t0,x0,y0,v0,xf,yf,g,nx,K,tau)

p0 = [x0;y0;v0;z(1);z(2);z(3)];
tf_G = z(4);
pTot0 = z(5:end);
pTot0 = reshape(pTot0,2*nx,K-1);

options = odeset('RelTol',1e-8);
option  = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);

E = [];
t = [];
p = [];

for k = 1:K
    if isequal(k,1)
    p0 = p0;
    else
    p0 = pTot0(:,k-1);
    end

    tauspan = [tau(k),tau(k+1)];
    [tout,pout] = ode113(@braOde,tauspan,p0,options,t0,tf_G,g)
    ptf = pout(end,:).';
    
    if k<K
	E = [E;ptf-pTot0(:,k)];
	%;Htf+1
    end
	t = [t; tout];
	p = [p; pout];
end
xtf     = ptf(1);
ytf     = ptf(2);
vtf     = ptf(3);
lamx_tf = ptf(4);
lamy_tf = ptf(5);
lamv_tf = ptf(6);
theta_G_tf = 0;
theta_tf = fsolve(@solveControl,theta_G_tf,option,lamx_tf,lamy_tf,lamv_tf,vtf,g);
Htf = lamx_tf*vtf*sin(theta_tf)+(lamy_tf*vtf+lamv_tf*g)*cos(theta_tf);
E = [E;ptf(1)-xf;ptf(2)-yf;lamv_tf;Htf+1];

end


%% state and costate function
function pdot = braOde(tspan,p,t0,tf,g)

theta_G = 0;
options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
theta = fsolve(@solveControl,theta_G,options,p(4),p(5),p(6),p(3),g);

x_dot = p(3)*sin(theta); % v*sind(theta)
y_dot = p(3)*cos(theta); % v*cosd(theta)
v_dot = g*cos(theta);
lamx_dot = 0;
lamy_dot = 0;
lamv_dot = -p(4)*sin(theta) - p(5)*cos(theta);

pdot = (tf-t0)/2*[ x_dot;y_dot;v_dot;lamx_dot;lamy_dot;lamv_dot];

end


%% H_theta function
function H_th = solveControl(theta,lamx,lamy,lamv,v,g)

H_th = lamx*v*cos(theta) - (lamy*v + lamv*g)*sin(theta); % function of theta

end