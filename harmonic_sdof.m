% harmonic_sdof.m
% harmonic response of the sdof mass-spring-damper
% plots Magnitude Ratio (dynamic amplification) and phase shift 
%         vs. excitation freq 
%
% last modified 11/02/18 CLee
%
close all
clear all

r = [0: 0.01: 2];             % r = omega_f/omega_n (forcing freq/natural freq)
zeta = [0.05 0.1 0.2 0.5 sqrt(2)/2  3.0];             % zeta = damping ratio

colors = {'b' 'g' 'r' 'k' 'c' 'm' 'y'}

for i = 1:length(zeta)
X = 1./ (sqrt( (1-r.^2).^2 + (2*zeta(i)*r).^2)) ;   % Displacement magnitude ratio
phi = atan2(2*zeta(i)*r, 1-r.^2);                   % Phase shift

figure(1)
plot(r,X, char(colors(i)))
hold on

figure(2)
plot(r,phi*180/pi,  char(colors(i)))                % phase shift in degrees
hold on 
end


figure(1) 
xlabel('Frequency Ratio')
ylabel('Magnitude Ratio')
title('Displacement Amplitude Ratio')
grid on
 

figure(2) 
xlabel('Frequency Ratio')
ylabel('Phase Shift (degrees)')
grid on
