%% main.m
function Brachistochrone_OCP
    close all;clear;clc;
    x0 = 0; y0 = 0; v0 = 0;t0 = 0;
    lamx =0;lamy = 0;lamv = 0;tf = 10;
    lam_guess = [lamx,lamy,lamv,tf];
    %P = [x0,y0,v0,lamx,lamy,lamv]';
    option = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);

    lam = fsolve(@Error_Brac,lam_guess,option,x0,y0,v0,t0);
    lam
end

%% ErrorFunction
function E = Error_Brac(lam_guess,x0,y0,v0,t0)
    xf = 2; 
    yf = 2;
    g = 10;
    tf = lam_guess(4);
    P0 = [x0;y0;v0;lam_guess(1);lam_guess(2);lam_guess(3)];
    
    option = odeset('RelTol',1e-8);
    [t,P] = ode113(@myOde_Brac,[t0 tf],P0,option,g);
    Ptf = P(end,:);
    xtf = Ptf(1); ytf = Ptf(2); vtf = Ptf(3);
    lamx_tf = Ptf(4); lamy_tf = Ptf(5); lamv_tf = Ptf(6);
    
    g = 10;
    options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
    theta_guess_tf = 0;
    theta_tf = fsolve(@solveControl,theta_guess_tf,options,lamx_tf,lamy_tf,lamv_tf,vtf,g);
    Htf=lamx_tf*vtf*sin(theta_tf)+(lamy_tf*vtf+lamv_tf*g)*cos(theta_tf);
    
    E = [xtf-xf; ytf-yf; lamv_tf; Htf+1]; 
    if E <= [1e-8; 1e-8; 1e-8; 1e-8]
        figure(1);
        plot(t,P(:,1),'r*-');
        title('x(t)');
        hold on
        grid on

        plot(t,P(:,2),'bo-');
        title('y(t)');
        hold on
        grid on

        plot(t,P(:,3),'ro-');
        title('v(t)');
        grid on
        hold on
        
        legend('x(t)','y(t)','v(t)');
    end

end

%% state and costate variables
function Pdot = myOde_Brac(t,P,g)
    options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
    
    theta_guess = 0;
    theta = fsolve(@solveControl,theta_guess,options,P(4),P(5),P(6),P(3),g);
    
    x_dot = P(3)*sin(theta); % v*sind(theta)
    y_dot = P(3)*cos(theta); % v*cosd(theta)
    v_dot = g*cos(theta);
    lamx_dot = 0;
    lamy_dot = 0;
    lamv_dot = -P(4)*sin(theta) - P(5)*cos(theta);

    Theta = [theta]
    
    
    figure(4);
    plot(t,theta,'bo-');
    title('theta(t)');
    hold on
    grid on;
        
    Pdot = [ x_dot;y_dot;v_dot;lamx_dot;lamy_dot;lamv_dot];
end

%% H_theta function
function H_th = solveControl(theta,lamx,lamy,lamv,v,g)

    H_th = lamx*v*cos(theta) - (lamy*v + lamv*g)*sin(theta);

end
