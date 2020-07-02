%% Adsorption Model

function Model
close all
clear
clc

%% Initialisation

t0 = 0;         % Initial Time s
tf = 10;       % Final time s
dt = .1;         % Time step
t = [t0:dt:tf]'; % Time vector

L = 50;         % Max Column length (m)
z = [0:.5:L]';   % Mesh generation
n = length(z);  % Size of mesh grid

c_const = 1;
%% ODE15S Solver

y_fun = @(x) heaviside(x-15).*(x-15) - heaviside(x-25).*(x-25) ...
    - heaviside(x-25).*(x-25) + heaviside(x-35).*(x-35);

y_fun = @(x) sin((pi/15)*x);
% Initial Conditions
% Initial condition for the entrance only (c(x=0,t=0) = 15, c(x!=0,t=0) = 0)
% c0 = zeros(n,1);
% c0(1) = 15;
% c0 = 25*ones(n,1); % Initial condition for all x postions (c(x,t=0) = 15)
c0 = y_fun(z);

figure
plot(z,y_fun(z))

y0 = [c0];
function dydt = deq(t, y)
    load parameters
    c = y(1:n);
    dcdt = zeros(n,1);
    dcdz = zeros(n,1);
    d2cdz2 = zeros(n,1);
    
    % Boundary condition for x0 at all times
    % c(x=0,t)
    dcdz(1) = 0;
    d2cdz2(1) = 0;
    
    for i=2:n-1
      dcdz(i) = (c(i)-c(i-1))/(z(i)-z(i-1));
      d2cdz2(i) = (c(i+1)-2*c(i)+c(i-1))/(z(i+1)-z(i))^2;
    end
    
    % Boundary condition for xn at all times
    dcdz(n) = 0;
    d2cdz2(n) = 0;

%     for i=1:n
%         dcdt(i) = -d2cdz2(i);
%     end
    % Differential Equation
    dcdt = -c_const*dcdz;
    
    dydt = [dcdt];
end

%% Numerical Solution

save parameters
options=odeset('RelTol',1e-9,'AbsTol',1e-9);

[t, Y] = ode15s(@deq,t,y0,options);
c_matrix = Y(:,1:n);
tot_matrix = Y;

% [X,Y] = meshgrid(z(1:1:end,1:1:end),t(1:1:end,1:1:end));
% Z = c_matrix(1:1:end,1:1:end);
% figure
% h = surf(X,Y,Z)
% set(h,'LineStyle','none')

[X,Y] = meshgrid(z,t);
Z = c_matrix;
figure
c = surf(X,Y,Z);
% set(c,'LineStyle','none')

title('Numerical solution computed with 20 mesh points')
xlabel('Distance x')
ylabel('Time t')


%% Analatical Solution
% Differential Equation dcdt = -c_const*dcdz;

[X,Y] = meshgrid(z,t);
Z = zeros(length(t),n);
for i = 1:1:length(t)
    Z(i,:) = y_fun(z-c_const*t(i));
end
figure
c = surf(X,Y,Z);
% set(c,'LineStyle','none')

title('Analatical solution')
xlabel('Distance x')
ylabel('Time t')



end