clear;

%% Simulation Settings
N = 3; %< number of timesteps to do reachability (index starts at 1)
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
xDim = @(k) memZono.genKeys(sprintf('x_k=%d',k),1:n);
uDim = @(k) {sprintf('u_k=%d',k)}; 

% Initial Conditions
X0 = zono(diag([1,2]),zeros(2,1));
X{1} = memZono(X0,xDim(1));  %<== lbl not needed since n=nG
Xall = X{1};

% Terminal Set
XF = zono(0.25*eye(2),ones(2.,1));
XF = memZono(XF,xDim(N)); %<== lbl not needed since n=nG

% Nominal Input Set
Unom = zono(1.5,0);

%% Reachability Calculation
switch 'method4' %'method1' 'method2' 'method3' 'method4'
    case 'method1' 
% Calculate and Save
for k = 1:N-1 % Time-evolution
    % Current Input
    U{k} = memZono(Unom,uDim(k));  %<== lbl not needed since n = nG
    % Step Update
    X{k+1} = X{k}.map(A,xDim(k),xDim(k+1)) ...
            + U{k}.map(B,uDim(k),xDim(k+1));
    % Save Data ($\starcross$)
    Xall = [Xall; U{k}; X{k+1}]; 
end
% Add Terminal Constraints
Xinter = Xall.and(XF, 'termCon');

    case 'method2' % recursively calcualte and projection
% map(), plus(), cartProd(), and()
for k = 1:N-1 % Time-Evolution
% Current Input
Uk = memZono(Unom,uDim(k)); %<== lbl not needed since n = nG
Xall = cartProd(Xall, Uk);
% Step Update
Xall = cartProd(Xall,...
    plus(Xall.map(A,xDim(k),xDim(k+1)),...
        Xall.map(B,uDim(k),xDim(k+1))));
end
% Add Terminal Constraints
Xinter = and(Xall,XF, 'termCon');
% Time-projections
for k = 1:N-1
    X{k} = Xall(xDim(k));
    U{k} = Xall(uDim(k));
end
X{N} = Xall(xDim(N));

for k=1:N; X{k}=Xall(xDim(k)); end
for k=1:N-1;U{k}=Xall(uDim(k));end

            
    case 'method3' % map(), plus(), cartProd(), and()
% time-Evolution
for k = 1:N-1
    % Current Input
    Xall = cartProd(Xall,memZono(Unom,uDim(k))); %<== lbl not needed since n = nG
    % Recursive Set Update
    Xall = cartProd(Xall,...
        plus(map(Xall,A,xDim(k),xDim(k+1)),...
            map(Xall,B,uDim(k),xDim(k+1))));
end
% Add Terminal Constraints
Xinter = and(Xall,XF,'termCon');
% Set Projections
X = arrayfun(@(k) Xall(xDim(k)),1:N,UniformOutput=false);
U = arrayfun(@(k) Xall(uDim(k)),1:N-1,UniformOutput=false);

    case 'method4' % Recursive Map and Save
% Time-Evolution
for k = 1:N-1
    % Current Input
    U{k} = memZono(Unom,uDim(k)); %<== lbl not needed since n = nG
    Xall = [Xall; U{k}];
    % Time-Update
    X{k+1} = Xall.map(A,xDim(k),xDim(k+1)) ...
            + Xall.map(B,uDim(k),xDim(k+1));
    Xall = [Xall; X{k+1}];
end
% Add Terminal Constraints
Xinter = and(Xall,XF,'termCon');
end

%% Plotting
fig = figure;

% State plots
subplot(1,2,1);
hold on;
plot(XF.Z(xDim(N)), 'g', 1);
drawnow;
for k = 1:N
    plot(Xinter.Z(xDim(k)), selectColor(k), 0.6);
    plot(X{k}.Z(xDim(k)), selectColor(k), 0.2);
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
uDims = [uDim(1),uDim(2)];
plot(Xall,uDims,'b',0.2);
plot(Xinter,uDims,'b',0.6);
% plot(Xall,[uDim(1),uDim(2)],'b',0.2);
% plot(Xinter,[uDim(1),uDim(2)],'b',0.6);
% plot(Z([U{1}; U{2}],[uDim(1),uDim(2)]),'b',0.2);
% plot(Xinter.Z([uDim(1),uDim(2)]),'b',0.6);
drawnow
hold off;

axis equal;
xlim([-2 2]);
ylim([-2 2]);
xlabel('$u(1)$','Interpreter','latex');
ylabel('$u(2)$','Interpreter','latex');



%%% Save Fig:
saveas(fig,'ex_reachability_basic.png')

%% local functions
function color = selectColor(i)
    colors = {'k','b','r'};
    color = colors{mod(i,length(colors))+1};
end


%% Old Code
% varDim = @(var,n,k) memZono.genKeys(sprintf('%s_k=%d',var,k),1:n);
% xDim = @(k) varDim('x',n,k);
% if m > 1, uDim = @(k) varDim('u',m,k); 
% else, uDim = @(k) {sprintf('u_k=%d',k)}; end

