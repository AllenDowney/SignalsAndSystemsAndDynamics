% laplace transform examples 
% using symbolic toolbox
%
% ENGR 2340 Dynamics 
% last modified 11/11/14 CLee

% define symbolic variables
syms s t


% Example 1
f1 = exp(-3*t)
F1 = laplace(f1)
G1    = exp(-3*s) / (s+2)^2
G11    = ilaplace(G1)
G11num = vpa(G11,4)
t1 = [0:  0.1: 10];
x1 = subs(G11num, t1);
figure (1)
plot(t1,x1)
title('Example 1')
xlabel ('time')
ylabel ('x(t)')


% Example 2
f2 = exp(-3*t)/((t+2)^2)
F2 = laplace(f2)
F2s = simplify(F2)
f2inv = ilaplace(F2s)
f2vpa = vpa (f2inv,4)
t2 = [0:  0.1: 1];
x2 = subs(f2vpa, t2);
figure (2)
plot(t2,x2)
title('Example 2')
xlabel ('time')
ylabel ('x(t)')

%Example 3
% example: f(t) = [ t*u(t)-2*(t-1)*u(t-1)+(t-2)*u(t-2) ]
C1 = sym ('Heaviside(t-1)')
C2 = sym('Heaviside(t-2)')

f3 = t - 2*(t-1)*C1 + (t-2)*C2
F3 = laplace(f3)
F3inv = ilaplace(F3)


%Example 4
N1 = [ 3 2 1];
D1 = [1 7 14 8];
[r,p,k]= residue(N1,D1)

%Example 5
N2 = [ 4 2];
D2 = [1 11 44 60];
[r,p,k]= residue(N2,D2)

%Example 6
N3 = [3 3 4];
D3 = [1 5 6];
[r,p,k]= residue(N3,D3)
