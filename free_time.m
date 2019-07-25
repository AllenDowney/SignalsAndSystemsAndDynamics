% free_time.m
%
% Step Response:
% single degree-of-freedom oscillator
% mass-spring-damper
%
% last modified 10/31/18 CLee
%
function free_time
clear all
close all
clear functions
%
% specify mass, spring constant, damping constant
% m =  1.0;            
% c =  0.0;
% k =  4.0;
% wn = sqrt(k/m);
% zeta = c/(2*wn*m);
% wn2 = k/m;

% or specify directly
wn = 4;
wn2 = wn*wn;
zeta =.1;
wd = wn*sqrt(1-zeta^2);
%
% roots of characteristic eqns (lambdas)
r1 = -zeta*wn + sqrt(zeta^2-1)*wn;
r2 = -zeta*wn - sqrt(zeta^2-1)*wn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time span
t_span = [0, 12];                     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state variables Z_1 = x, Z_2 = x_dot,
x0 =  2;             % initial displacement
v0 =  0;             % initial velocity
Z_0 = [x0, v0];      
% 
reltol = 1.0e-8;
options= odeset('RelTol', reltol);
[t, zout] = ode113(@sdof_fun, t_span, Z_0, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical solutions
% check zeta values
% undamped step response with IC's = 0 ***
    if  zeta == 0                                                 %undamped
xanalytic = v0/wn*sin(wn*t) + (x0)*cos(wn*t);
    elseif zeta > 1                                            %overdamped
xanalytic = (x0*r2-v0)/(r2-r1)*exp(r1*t) + (x0*r1-v0)/(r1-r2)*exp(r2*t);    
    elseif zeta == 1                                     %critically damped
xanalytic = exp(-wn*t).*( (v0+x0*wn)*t + x0 );   
    elseif  zeta < 1 & zeta > 0                                %underdamped
xanalytic = exp(-zeta*wn*t).*( x0*cos(wd*t) +  (v0+x0*zeta*wn)/wd*sin(wd*t));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure  
subplot(2,1,1)
plot( t, zout(:,1) )    
hold
plot(t,xanalytic,'r+')
xlabel('Time')
ylabel('Displacement')
title('SDOF Free Response with ICs')
legend('numerical','analytic')
%
subplot(2,1,2)
plot( t, zout(:,2) )      
xlabel('Time ')
ylabel('Velocity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOM's in state-space form
function dzdt = sdof_fun(T, ZZ)
% sdof oscillator equation in state space form
% second order oscillator as first order, state space form
zdot1 = ZZ(2);
zdot2 = -wn2*ZZ(1) - 2*zeta*wn*ZZ(2); 
% 
dzdt = [zdot1;zdot2];
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

