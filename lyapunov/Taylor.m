syms x1 x2 x3 x4
% x xdot theta thetadot

m = 1;
M = 5;
L = 2;
g = -10;
d = 1;
b = 1;

Sx = sin(x3);
Cx = cos(x3);
D = m*L*L*(M+m*(1-Cx^2));

%%  Design LQR controller

A = [0 1 0 0;
    0 -d/M b*m*g/M 0;
    0 0 0 1;
    0 -b*d/(M*L) -b*(m+M)*g/(M*L) 0];
B = [0; 1/M; 0; b*1/(M*L)];


Q = [10 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 1];
R = 1;

[K,S,e] = lqr(A,B,Q,R);

x = [x1;x2;x3;x4];
wr = [0; 0; pi; 0];      % reference position
u=-K*(x - wr);       % control law

f1 = x2;
f2 = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x4^2*Sx - d*x2)) + m*L*L*(1/D)*u;
f3 = x4;
f4 = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*(x4^2)*Sx - d*x2)) - m*L*Cx*(1/D)*u;

T1 = taylor(f1, [x1, x2, x3, x4], [0,0,pi,0],'Order',4);
T2 = taylor(f2, [x1, x2, x3, x4], [0,0,pi,0],'Order',4);
T3 = taylor(f3, [x1, x2, x3, x4], [0,0,pi,0],'Order',4);
T4 = taylor(f4, [x1, x2, x3, x4], [0,0,pi,0],'Order',4);

P = (x-wr)*(x-wr).';
V = dot(S, P)