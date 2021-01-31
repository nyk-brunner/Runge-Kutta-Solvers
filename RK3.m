function [x, y] = RK3(ode, a, b, h, Y, type, opt)
% RK3 uses the third order Runge-Kutta (roon-geh koot-tuh) method to solve 
% first-order ordinary differential equations.
%
% Inputs:
% First-order ordinary differential equation, ode
%   Left x-bound, a
%   Right x-bound, b
%   Step size, h
%   Value of solution at y (initial value), Y
%   What method type you'd like to use: (default: classical)
%       "Classical" or "Classic" for the classical R-K method
%       "Nystrom" for Nystrom's method
%       "Optimal" for the "nearly optimal" method
%       "Heun" for Heun's Third method
%
% Outputs:
%   Vector containing "independent" variable solution, x
%   Vector containing "dependent" variable solution, y
%
% Syntax:
%   [x, y] = RK3(ode, a, b, h, Y, type, opt)

if nargin == 5 % Default settings (Classical method, no plotting)
    type = 'classical';
    opt = 'no';
end

if strcmpi(type, 'plot') % Use default setting (Classical) and plot
    type = 'huen'; 
        fprintf('Default setting used (Classical) \n')
    opt = 'plot';
elseif nargin == 6
    opt = 'no';
end

% Initialize:
N = (b - a)/h;
x = zeros(N+1, 1);
y = zeros(N+1, 1);

% Apply initial condition:
x(1) = a;
y(1) = Y;

if strcmpi(type, 'classical') || strcmpi(type, 'classic')
    c1 = 1/6;
    c2 = 4/6;
    c3 = 1/6;
    a2 = 1/2;
    a3 = 1;
    b21 = 1/2;
    b31 = -1;
    b32 = 2;
elseif strcmpi(type, 'nystrom')
    c1 = 2/8;
    c2 = 3/8;
    c3 = 3/8;
    a2 = 2/3;
    a3 = 2/3;
    b21 = 2/3;
    b31 = 0;
    b32 = 2/3;
elseif strcmpi(type, 'optimal')
    c1 = 2/9;
    c2 = 3/9;
    c3 = 4/9;
    a2 = 1/2;
    a3 = 3/4;
    b21 = 1/2;
    b31 = 0;
    b32 = 3/4;
elseif strcmpi(type, 'Heun')
    c1 = 1/4;
    c2 = 0;
    c3 = 3/4;
    a2 = 1/3;
    a3 = 2/3;
    b21 = 1/3;
    b31 = 0;
    b32 = 2/3;
end

for i = 1:N
    x(i+1) = x(i) + h;
    K1 = ode(x(i),y(i));
    K2 = ode(x(i) + a2*h, y(i) + b21*K1*h);
    K3 = ode(x(i) + a3*h, y(i) + b31*K1*h + b32*K2*h);
    y(i+1) = y(i) + (c1*K1 + c2*K2 + c3*K3)*h;
end

if strcmpi(opt, 'plot')
    figure(1)
        plot(x,y,'k-','linewidth',1.5)
        
        grid on
        xlabel('X')
        ylabel('Y')
        title('ODE')
end