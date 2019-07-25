% laplace_ode.m 
% solves a single dof, 2nd-order ODE using Laplace transforms/Inverse LT's
% IC's are accounted for and can be non-zero
%
% last modified 11/8/18 CLee
%
clear all, close all

syms s t X  x x(t)  sX 

% system parameters  
% zeta  = 0.01;
% wn    = 10;
% % forcing parameters
% wf    = 12;
% A     = 1;

% IC's
x0    = 0;
v0    = 0;

% define ODE
% LHS = 2nd order ODE, RHS = forcing function 
% edit RHS to create desired forcing
%
% ode_eqn = '1.0e3*D(D(x))(t) + 1.0e2*D(x)(t) + 1.0e4*x(t)= 0.5*heaviside(t)'
% ode_eqn = '1.0e3*D(D(x))(t) + 1.0e2*D(x)(t) + 1.0e4*x(t)= 1.0e4*0.5*heaviside(t)+ 0.5*1.0e2*dirac(t)'

ode_eqn = diff( diff(x,t),t ) + 2* diff(x,t)+ x == sin(2*t)

% take the Laplace transform of both sides
s_ode =  laplace(ode_eqn,t,s)

% 
% % substitute in numerical values for system and forcing parameters and IC's
% % s_eqn = subs(s_ode,{'laplace(x(t), t, s)', 'x(0)', 'D(x)(0), zeta, wn^2, wn, wf, A'},...
% %              {X, x0, v0, zeta, wn^2, wn, wf, A} );
s_eqn = subs(s_ode,{laplace(x(t), t, s), subs(diff(x(t), t), t, 0), x(0)},...
             {X, v0, x0 } )

% solve equations for X(s)
X = solve(s_eqn,X)
% 
% % take inverse LT to get x(t) 
x = ilaplace(X)
% 
%  
% %plotting
t = linspace(0, 50, 500);
xplot = subs(x,t);   
plot(t,xplot, 'r')
xlabel('time')
ylabel('x(t)')



