clear;

%% Simulation Settings
N = 3; %< number of timesteps to do reachability
dt = 0.25;

%% System Definition
% DT-LTI System
A = [   1, 0.25;
     -0.5, 0.75];
B = [   0.025;
        0.25];
n = size(A,1); m = size(B,2);

%% Reachability Setup
% Setup Labels
xDims = @(k) memZono.genKeys(sprintf('%s_k=%d',var,k),1:n);
uDims = @(k) {sprintf('u_k=%d',k)}; 

% Initial Conditions
X_0 = zono(diag([1,2]),zeros(2,1));
X_{1} = memZono(X_0,xDims(1));
X_all = X_{1};

% Terminal Set
X_F = zono(0.25*eye(2),ones(2.,1));
X_F = memZono(X_F,xDims(N));

% Nominal Input Set
U_nom = zono(1.5,0);

%% Reachability Calculation
switch 'method4' %'method1' 'method2' 'method3' 'method4'
    case 'method1' 
% Calculate and Save
for k = 1:N-1 % Time-evolution
    % Current Input
    U_{k} = memZono(U_nom,uDims(k));
    % Step Update
    X_{k+1} = X_{k}.map(A,xDims(k),xDims(k+1)) ...
            + U_{k}.map(B,uDims(k),xDims(k+1));
    % Save Data ($\starcross$)
    X_all = [X_all; U_{k}; X_{k+1}]; 
end
% Add Terminal Constraints
X_inter = X_all.and(X_F, 'terminal_cons');

    case 'method2' % recursively calcualte and projection
% map(), plus(), cartProd(), and()
for k = 1:N-1 % Time-Evolution
% Current Input
U_k = memZono(U_nom,uDims(k));
X_all = cartProd(X_all, U_k);
% Step Update
X_all = cartProd(X_all,...
    plus(X_all.map(A,xDims(k),xDims(k+1)),...
        X_all.map(B,uDims(k),xDims(k+1))));
end
% Add Terminal Constraints
X_inter = and(X_all,X_F, 'terminal_cons');
% Time-projections
for k = 1:N-1
    X_{k} = X_all(xDims(k));
    U_{k} = X_all(uDims(k));
end
X_{N} = X_all(xDims(N));

for k=1:N; X_{k}=X_all(xDims(k)); end
for k=1:N-1;U_{k}=X_all(uDims(k));end

            
    case 'method3' % map(), plus(), cartProd(), and()
% time-Evolution
for k = 1:N-1
    % Current Input
    X_all = cartProd(X_all,memZono(U_nom,uDims(k)));
    % Recursive Set Update
    X_all = cartProd(X_all,...
        plus(map(X_all,A,xDims(k),xDims(k+1)),...
            map(X_all,B,uDims(k),xDims(k+1))));
end
% Add Terminal Constraints
X_inter = and(X_all,X_F,'terminal_cons');
% Set Projections
X_ = arrayfun(@(k) X_all(xDims(k)),1:N,UniformOutput=false);
U_ = arrayfun(@(k) X_all(uDims(k)),1:N-1,UniformOutput=false);

    case 'method4' % Recursive Map and Save
% Time-Evolution
for k = 1:N-1
    % Current Input
    U_{k} = memZono(U_nom,uDims(k));
    X_all = [X_all; U_{k}];
    % Time-Update
    X_{k+1} = X_all.map(A,xDims(k),xDims(k+1)) ...
            + X_all.map(B,uDims(k),xDims(k+1));
    X_all = [X_all; X_{k+1}];
end
% Add Terminal Constraints
X_inter = and(X_all,X_F,'terminal_cons');
end

%% Plotting
fig = figure;

% State plots
subplot(1,2,1);
hold on;
plot(X_F, 'g', 1);
drawnow;
for k = 1:N
    plot(X_inter(xDims(k)), selectColor(k), 0.6);
    plot(X_{k}, selectColor(k), 0.2);
    drawnow;
end
hold off;

axis equal;
xlim([-3 3]);
ylim([-3 3]);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');


% Input Plots
subplot(1,2,2);
hold on;
plot([U_{1}; U_{2}],'b',0.2);
plot(X_inter([uDims(1),uDims(2)]),'b',0.6);
drawnow
hold off;

axis equal;
xlim([-2 2]);
ylim([-2 2]);
xlabel('$u(1)$','Interpreter','latex');
ylabel('$u(2)$','Interpreter','latex');

%% local functions
function color = selectColor(i)
    colors = {'k','b','r'};
    color = colors{mod(i,length(colors))+1};
end


%% Old Code
% varDims = @(var,n,k) memZono.genKeys(sprintf('%s_k=%d',var,k),1:n);
% xDims = @(k) varDims('x',n,k);
% if m > 1, uDims = @(k) varDims('u',m,k); 
% else, uDims = @(k) {sprintf('u_k=%d',k)}; end

