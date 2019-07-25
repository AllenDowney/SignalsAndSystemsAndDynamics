% step_time_parameters.m
%
% Free Response for various zeta and omega_n
% single degree-of-freedom oscillator
% mass-spring-damper
%
% last modified 10/30/14 CLee
%
function step_time_parameters
clear all
close all
clear functions
%
% specify mass, spring constant, damping constant

% define step magnitude
A = 1;
% define time space
t_span = [0, 50];

% specify natural frequency and damping ratio arrays
wnarray = [0.25   0.5   1   2   10];
wn2array=wnarray.*wnarray;
zetaarray = [0.01    0.2  0.5  1.0   2.0];           

colors = {'b' 'g' 'r' 'k' 'c' 'm' 'y'};

% state variables Z_1 = x, Z_2 = x_dot,
x0 =  .0;   % initial displacement
v0 =  -0.0;   % initial velocity
Z_0 = [x0, v0];     %specify initial conditions (rads, rads/sec)
% 
reltol = 1.0e-8;
options= odeset('RelTol', reltol);

wn = wnarray(2);
for i = 1:length(zetaarray)

zeta = zetaarray(i);
[t, zout] = ode45(@sdof_fun, t_span, Z_0, options);

figure(1)
 plot(t,zout(:,1)*wn*wn, char(colors(i) ))
xlabel('Time') 
ylabel('Normalized Displacement, x(t)')
hold on
title('fixed nat. freq., variable damping ratio. with nonzero ICs')
end
%******************************
zeta = zetaarray(2);
for i = 1:length(wnarray)
wn = wnarray(i);
[t, zout] = ode45(@sdof_fun, t_span, Z_0, options);

figure(2)
plot(t,zout(:,1)*wn*wn, char(colors(i) ))
xlabel('Time') 
ylabel('Normalized Displacement, x(t)')
hold on
title('fixed damp ratio, variable nat. freq. with nonzero ICs')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOM's in state-space form
function dzdt = sdof_fun(T, ZZ)
% sdof oscillator equation in state space form
% second order oscillator as first order, state space form
dz1dt = ZZ(2);
dz2dt = -wn*wn*ZZ(1) - 2*zeta*wn*ZZ(2) + A; 
% 
dzdt = [dz1dt;dz2dt];
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

