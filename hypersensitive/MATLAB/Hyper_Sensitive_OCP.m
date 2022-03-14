%% main.m
function Hyper_Sensitive_OCP

    close all;clear;clc;
    x0 = 1;
    xf = 1;
    t0 = 0;
    tf = 10;
    lam_guess = 0;

    option = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
    lam0 = fsolve(@myError,lam_guess,option,x0,t0,xf,tf)

end

%% Error_function
function E= myError(lam0,x0,t0,xf,tf)

    P0 = [x0;lam0];
    option = odeset('RelTol',1e-8);
    tspan = [t0 tf];
    [t,P] = ode113(@myOde_Hyper,tspan,P0,option)
    Ptf = P(end,:)
    xtf = Ptf(1)
    E = xtf-xf;

    if E<=1e-8
        figure(1);
        plot(t,P(:,1),'r*-');
        title('x(t)');
        grid on

        figure(2);
        plot(t,-P(:,2),'bo-');
        title('u(t)');
        grid on
    end
    
end

%% state and costate function
function P_dot = myOde_Hyper(t,x)
    
    P_dot = [-1,-1;-1,1]*[x(1);x(2)];

end
