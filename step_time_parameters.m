% step_time_parameters.m
%
% Free Response for various zeta and omega_n
% single degree-of-freedom oscillator
% mass-spring-damper
%
% last modified 10/30/14 CLee
%
function step_time_parameters
%
% specify mass, spring constant, damping constant

% define step magnitude
A = 1;

% define time span
t_span = [0, 50];

% specify natural frequency and damping ratio arrays
wnarray = [0.25   0.5   1   2   10];
zetaarray = [0.01    0.2  0.5  1.0   2.0];           

% state variables Z_1 = x, Z_2 = x_dot
x0 =  0.0;    % initial displacement
v0 =  0.0;    % initial velocity
Z_0 = [x0, v0];     % initial conditions (rads, rads/sec)

% set graphics defaults
set(groot, 'DefaultLineLineWidth', 2);

% fix wn and sweep zeta
figure(1)
hold on

wn = wnarray(2);
for i = 1:length(zetaarray)
    zeta = zetaarray(i);
    [T, zout] = ode45(@sdof_fun, t_span, Z_0);

    label = sprintf('zeta=%0.2f', zeta);
    plot(T, zout(:,1)*wn*wn, 'DisplayName', label)
    xlabel('Time') 
    ylabel('Normalized Displacement, x(t)')
    title('fixed nat. freq., variable damping ratio. with nonzero ICs')
    legend()
end

% fix zeta and sweep wn
figure(2)
hold on

zeta = zetaarray(2);
for i = 1:length(wnarray)
    wn = wnarray(i);
    [T, zout] = ode45(@sdof_fun, t_span, Z_0);

    label = sprintf('wn=%0.2f', wn);
    plot(T, zout(:,1)*wn*wn, 'DisplayName', label)
    xlabel('Time') 
    ylabel('Normalized Displacement, x(t)')
    title('fixed damp ratio, variable nat. freq. with nonzero ICs')
    legend()
end


function dzdt = sdof_fun(t, Z)
    % sdof second order oscillator as first order, state space form
    dz1dt = Z(2);
    dz2dt = -wn*wn*Z(1) - 2*zeta*wn*Z(2) + A; 
    
    dzdt = [dz1dt; dz2dt];
end

end

