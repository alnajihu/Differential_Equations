%% Adsorption Model

function Heat_Equation
close all
clear
clc

%% Initialisation

L = 1;
z = linspace(0,L,20);
t = [linspace(0,0.05,20), linspace(0.5,5,10)];
n = length(z);          % Size of mesh grid

c_const = 1;

%% Heat Solver

% Initial condition for all x postions (c(x,t=0) = 0.5)
c0 = 0.5*ones(n,1);
c0(1) = 0;
c0(n) = 1;

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
    
    % Boundary condition for xn at all times
    % c(x=n,t)
    dcdz(n) = 0;
    d2cdz2(n) = 0;
    
    for i=2:n-1
      dcdz(i) = (c(i)-c(i-1))/(z(i)-z(i-1));
      d2cdz2(i) = (c(i+1)-2*c(i)+c(i-1))/(z(i+1)-z(i))^2;
    end

    % Differential Equation
    dcdt = c_const*d2cdz2;
    
    dydt = [dcdt];
end

%% Numerical Solution

save parameters
options=odeset('RelTol',1e-9,'AbsTol',1e-9);

[t, Y] = ode15s(@deq,t,y0,options);
c_matrix = Y(:,1:n);
tot_matrix = Y;

[X,Y] = meshgrid(z,t);
Z = c_matrix;

%%
% figure
% colormap hot
% c = surf(X,Y,Z);
% % set(c,'LineStyle','none')
% 
% title('Numerical solution computed with 20 mesh points')
% xlabel('Distance x')
% ylabel('Time t')



%%
figure
colormap hot
imagesc(z,t,Z)
colorbar
xlabel('Distance x','interpreter','latex')
ylabel('Time t','interpreter','latex')
title('Heat Equation for $0 \le x \le 1$ and $0 \le t \le 5$','interpreter','latex')


%% Analatical Solution
% % Differential Equation dcdt = -c_const*dcdz;
% 
% [X,Y] = meshgrid(z,t);
% Z = zeros(length(t),n);
% for i = 1:1:length(t)
%     Z(i,:) = y_fun(z-c_const*t(i));
% end
% figure
% c = surf(X,Y,Z);
% % set(c,'LineStyle','none')
% 
% title('Analatical solution')
% xlabel('Distance x')
% ylabel('Time t')



end