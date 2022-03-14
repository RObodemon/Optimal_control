%% main function
function hyper_Sensitive_OCP_main

close all;clear;clc

nx = 1; % # of states
x0 = 1;
xf = 1;
t0 = 0;
tf = 50; %10, 20, 30, 40, 50
K  = 3;     % 2, 4, 8, 16, 20

tau = linspace(-1,+1,K+1);
lambda0_G = 0;
P0_G      = zeros(2*nx,K-1);
z_G       = [lambda0_G; P0_G(:)];       % into column

options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
z = fsolve(@hyperSensitiveError,z_G,options,x0,t0,xf,tf,nx,K,tau)
[E,t,p] = hyperSensitiveError(z,x0,t0,xf,tf,nx,K,tau);

% plot the trajectory
figure(1);
plot(t,p(:,1),'r*-');
title('x(t) / tf=50,k=3');
grid on;

% plot the control
figure(2);
plot(t,-p(:,2),'bo-');
title('u(t) / tf=50,k=3');
grid on;

end


%% Error function

function [E,t,p] = hyperSensitiveError(z,x0,t0,xf,tf,nx,K,tau)  %variables and constraints

lambda0 = z(1);
pTot0   = z(2:end)
pTot0   = reshape(pTot0,2*nx,K-1);

options = odeset('RelTol',1e-8);

E = [];
t = [];
p = [];

for k =1:K
    if isequal(k,1)
        p0 = [x0;lambda0];
    else
        p0 = pTot0(:,k-1); % k-1 roll and all column
    end

    tauspan     = [tau(k),tau(k+1)];
    [tout,pout] = ode113(@hyperSensitiveODE,tauspan,p0,options,t0,tf);
    ptf         = pout(end,:).'; %final result

    if k<K
        E = [E;ptf(1)-pTot0(1,k)];
    end
    t = [t; tout];
    p = [p; pout];
end
E = [E; ptf(1)-xf];

end


%% state and costate function

function pdot = hyperSensitiveODE(tspan,p,t0,tf)

A    = [-1 -1; -1 1];
pdot = (tf-t0)/2*A*p; % scale the dynamic

end