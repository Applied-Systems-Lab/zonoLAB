% Programmer Name: Justin Chen
% Problem Definition: Given some initial state utilize zonotopes to evolve
% the state into the final state, a targeted state. Then show the possible
% input values that would allow the initial state to reach the final state.

clear;
%clc; close all;

%% Simulation Settings
T = 3;
dt = 1.5;
N = round(T/dt)+1;
time = linspace(0,T,N);

%% System Definition
% Continuous time dynamics
A = [   0, 1;
        0, 0];

B = [   0;
        1];

x0 = [0; 0];
%% Discretization
% Discrete Time dynamics
F = expm(A*dt);
[n,m] = size(B);

tau = linspace(0,dt,1001);
expmAtauB = zeros(n,m, length(tau));

for i = 1: length(tau)
    expmAtauB(:,:,i) = expm(A*tau(i))*B;
end

G = trapz(tau, expmAtauB, 3);   % numerical integration across 3 dimensions

%%
% Input Parameters
Gu = 1;     % Range
Cu = 1;     % Offset

% Starting Position
Gx = [1 0;      % Col Sum length
      0 1];
Cx = x0;    % Center

% Final Position
Gf = [  1 0;    % Col Sum length
        0 1];
Cf = [10;5];   % Center

%% Reachability Setup
X_0 = zono(Gx,Cx);      % Make the initial state sized by first parameter and located at x0
X_F = zono(Gf,Cf);              % Make the final state sized by Gx and located at Cx
U_nom = zono(Gu,Cu);         % Zono based on G = u and C = offset

%% Reachability Calculation
X_{1} = memZono(X_0,'x_1');
X_all = X_{1};

X_F = memZono(X_F, sprintf('x_%d',N));

for k = 1:N-1
    U_{k} = memZono(U_nom, sprintf('u_%d',k));

    newDims = {sprintf('x_%d_1',k+1),sprintf('x_%d_2',k+1)};
    X_{k+1} = X_{k}.transform(U_{k}.transform([],G,{},newDims),F,{},newDims);

    % Stacking
    X_all = X_all.merge(U_{k});
    X_all = X_all.merge(X_{k+1});
end

X_inter = X_all.merge(X_F,'terminal_cons'); % Intersect Common Dimension

%% Plotting
fig = figure;

% State Plot
subplot(1,2,1);
hold on;
plot(X_F, 'all', 'k', 0.1);
drawnow;

for k = 1:N
    plot(X_{k}, 'all', selectColor(k), 0.2)
    plot(X_inter, {sprintf('x_%d_1',k),sprintf('x_%d_2',k)}, selectColor(k), 0.6);
end
hold off;

xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');

%% Input Plots
subplot(1,2,2);
hold on;
plot([U_{1}; U_{2}],'all','b',0.2);
plot(X_inter,{'u_1','u_2'},'b',0.6);

xlabel('$u_1$','Interpreter','latex');
ylabel('$u_2$','Interpreter','latex');

axis equal
drawnow

%% Paths from intial state
x_const = [0.5;0.5];

G_C = zono(zeros(n),x_const);

G_{1} = memZono(G_C,'x_1');
G_all = G_{1};

for k = 1:N-1
    newDims = {sprintf('x_%d_1',k+1),sprintf('x_%d_2',k+1)};
    G_{k+1} = G_{k}.transform(U_{k}.transform([],G,{},newDims),F,{},newDims);

    % Stacking
    G_all = G_all.merge(U_{k});
    G_all = G_all.merge(G_{k+1});
end

G_inter = G_all.merge(X_F,'terminal_cons'); % Intersect Common Dimension

subplot(1,2,1);
hold on;

for k = 1:N
    plot(G_inter, {sprintf('x_%d_1',k),sprintf('x_%d_2',k)}, selectColor(k+1), 1);
%   plot(G_{k}, 'all', selectColor(k+1), 0.2)
end

hold off;

%% Input Plots for path
subplot(1,2,2);
hold on;
plot(G_inter,{'u_1','u_2'},'k',1);
drawnow
hold off

%% Testing

subplot(1,2,1);
hold on;
x_const = [0.5;0.5];

% u is a list of inputs, the graph u_1 and u_2 is the list of inputs at
% each iteration.

u_guess = [5/3,2];    % At boundary
% u_guess = [1.8,2];    % Slightly outside boundary
% 
for k = 2:N
    plot(x_const(1),x_const(2),'kx');
    x_const = F*x_const + G*u_guess(k-1);
end

    plot(x_const(1),x_const(2),'kx');

subplot(1,2,2);
hold on;
plot(u_guess(1),u_guess(2),'wx');
hold off

% % 100 samples
% for count = 1:1000
%     subplot(1,2,1);
%     hold on;
%     x_const = [0.5;0.5];
%     plot(x_const(1),x_const(2),'ko');
% 
%     u_guess = (2*Gu).*rand(1,n);    % Slightly outside boundary
%     
%     for k = 2:N
%         x_const = F*x_const + G*u_guess;
%         plot(x_const(1),x_const(2),'k.');
%     end
% 
%     subplot(1,2,2);
%     hold on;
%     plot(u_guess(1),u_guess(2),'wx');
% end
% hold off

%% local functions
function color = selectColor(i)
    colors = {'k','r','b','c','m'};
    color = colors{mod(i,length(colors))+1};
end