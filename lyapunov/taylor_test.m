% Get the desired symbolic expression
% syms x1 x2

% g = 9.81;
% m = 1.;
% l = 0.5;
% d = 0.1;

% convert symbolic expression to function handle
% gFH = matlabFunction(g);

% define dynamics
syms x1 x2 g m l d k1 k2;

u = -(x1-pi)*k1 - k2*x2;
f1 = x2;
f2 = (u-d*x2-sin(x1)*m*g*l)/(m*l*l);

% 1 st
T = jacobian([f1], [x1, x2])

T20 = hessian(f1, [x1, x2])
T21 = hessian(f2, [x1, x2])




diff((g*sin(x1))/l, x1)






