function [x, y] = RK2(ode, a, b, h, Y, type, opt)
% RK2 uses various second-order Runge-Kutta(roon-geh - koot-tuh) methoda
% to solve a first order ODE. Methods included are the modified Euler
% (in the form of a second order R-K), the midpoint method, and Huen's
% method.
%
% Inputs:
%   First-order ordinary differential equation, ode
%   Left x-bound, a
%   Right x-bound, b
%   Step size, h
%   Value of solution at y (initial value), Y
%   What method you would like to use: (default: Huen)
%       "Euler" for modified Euler method (R-K converted)
%       "Mid" for midpoint method
%       "Huen" for Huen's method
%   Plotting option, opt
%       If you want to the function to plot, set opt = plot
%
% Outputs:
%   Vector containing "independent" variable solution, x
%   Vector containing "dependent" variable solution, y
%
% Syntax:
%   [x, y] = RK2(ode, a, b, h, Y, type, opt)

if nargin == 5
    type = 'huen';
    opt = 'no';
end

if strcmpi(type, 'plot') % Use default setting (Huen's)
    type = 'huen'; 
        fprintf('Default setting used (Huen''s) \n')
    opt = 'plot';
elseif nargin == 6
    opt = 'no';
end

if strcmpi(type,'euler') % Type selected is the modified Euler (R-K form)
    c1 = 1/2;
    c2 = 1/2;
    a2 = 1;
    b21 = 1;
elseif strcmpi(type, 'mid') % Type selected is the midpoint method
    c1 = 0;
    c2 = 1;
    a2 = 1/2;
    b21 = 1/2;
elseif strcmpi(type, 'huen') % Type selected is Huen's method
    c1 = 1/4;
    c2 = 3/4;
    a2 = 2/3;
    b21 = 2/3;
end

% Initialize:
N = (b - a)/h; % Number of steps
x = zeros(N+1,1);
y = zeros(N+1,1);

% Apply conditions:
x(1) = a;
y(1) = Y;

for i = 1:N
    x(i+1) = x(i) + h; % Update x
    K1 = ode(x(i), y(i));
    K2 = ode(x(i) + a2*h, y(i) + b21*K1*h);
    y(i+1) = y(i) + (c1*K1 + c2*K2)*h;
end

% Plot:
if strcmpi(opt, 'plot')
    figure(1)
        plot(x,y,'b-','linewidth',1.5)
        
        grid on
        xlabel('X')
        ylabel('Y')
        title('ODE')
end
