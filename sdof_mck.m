% sdof_mck.m
%
% Single degree-of-freedom oscillator
% Mass-Spring-Damper
%
% last modified 1/28/19 CLee
%
function sdof_mck
clear all
close all
clear functions
%
% specify mass, spring constant, damping constant
m =  1.0;             % units ?
c =  0.0;
k =  4.0;
 
wn = sqrt(k/m);
zeta = c/(2*wn*m);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define forcing function: forcef
% can be general function, harmonic, or free

% Define simulation parameters
t_span = [0:0.1:50];          % max time span for simulation 

% step response steady-state value
magstep = 20;
ss = magstep/wn^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state variables Z_1 = x, Z_2 = x_dot,
x0 = 0.0;   % initial displacement
v0 = 0.0;   % initial velocity
Z_0 = [x0, v0];     %specify initial conditions (rads, rads/sec)
% 
reltol = 1.0e-10;
options= odeset('RelTol', reltol);
[t, zout] = ode113(@sdof_fun, t_span, Z_0, options);

% analytical soln: undamped step response with IC's = 0 ***
xundamp = v0/wn*sin(wn*t) + (x0-magstep/wn^2)*cos(wn*t) + magstep/wn^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure  % angular displacement, velocity, and acceleration
subplot(2,1,1)
plot( t, zout(:,1) )    
hold
plot(t,xundamp,'r+')
xlabel('Time (sec)')
ylabel('Displacement (length)')
title('SDOF')

subplot(2,1,2)
plot( t, zout(:,2) )      
xlabel('Time (sec)')
ylabel('Velocity (length/sec)')
title('SDOF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOM's in state-space form
function dzdt = sdof_fun(T, Z)
% sdof oscillator equation in state space form

% Various forcing functions on RHS
%
% for Step function with time delay
% general functions: step with time shift
% if T < 0
%   forcef = 0;
%    %forcef = 25/10*t;
%    %forcef = -100/20*t+100;
% else
%    forcef = magstep;
%   %forcef = 25;
%   %forcef=0;
% end

% for Harmonic excitation
a = 1;
forcef= a*sin(2*T);

% ffor Free vibration
% forcef= 0;

% second order oscillator in first order, state space form
dz1dt = Z(2);
dz2dt = -k*Z(1)/m - c/m*Z(2) + forcef/m; 

dzdt = [dz1dt;dz2dt];
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

