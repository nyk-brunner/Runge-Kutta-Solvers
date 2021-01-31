function [x, y] = RK4(ode, a, b, h, Y, opt)
% RK4 uses the fourth-order Runge-Kutta (roon-geh koot-tuh) method of
% solving first-order ODE's.
%
% Inputs:
% First-order ordinary differential equation, ode
%   Left x-bound, a
%   Right x-bound, b
%   Step size, h
%   Value of solution at y (initial value), Y
%
% Outputs:
%   Vector containing "independent" variable solution, x
%   Vector containing "dependent" variable solution, y
%
% Syntax:
%   [x, y] = RK3(ode, a, b, h, Y, opt)

if nargin == 5
    opt = 'no';
end

% Initialize:
N = (b - a)/h;
x = zeros(N+1, 1);
y = zeros(N+1, 1);

% Apply initial value:
x(1) = a;
y(1) = Y;

% Constants:
c1 = 1/6;
c2 = 2/6;
c3 = 2/6;
c4 = 1/6;
a2 = 1/2;
a3 = 1/2;
a4 = 1;
b21 = 1/2;
b31 = 0;
b32 = 1/2;
b41 = 0;
b42 = 0;
b43 = 1;

for i = 1:N
    x(i+1) = x(i) + h;
    K1 = ode(x(i), y(i));
    K2 = ode(x(i) + a2*h, y(i) + b21*K1*h);
    K3 = ode(x(i) + a3*h, y(i) + b31*K1*h + b32*K2*h);
    K4 = ode(x(i) + a4*h, y(i) + b41*K1*h + b42*K2*h + b43*K3*h);
    y(i+1) = y(i) + (c1*K1 + c2*K2 + c3*K3 + c4*K4)*h;
end

if strcmpi(opt, 'plot')
    figure(1)
        plot(x,y,'k-','linewidth',1.5)
        
        grid on
        xlabel('X')
        ylabel('Y')
        title('ODE')
end
