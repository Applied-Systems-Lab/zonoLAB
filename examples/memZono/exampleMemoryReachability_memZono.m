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
X_0 = zono(diag([1,2]),zeros(2,1));
X_F = zono(0.25*eye(2),ones(2.,1));
U_nom = zono(1.5,0);

%% Reachability Calculation
% Setup Labels
varDims = @(var,n,k) memZono.genKeys(sprintf('%s_k=%d',var,k),1:n);
xDims = @(k) varDims('x',n,k);
if m > 1, uDims = @(k) varDims('u',m,k); 
else, uDims = @(k) {sprintf('u_k=%d',k)}; end

% Initial Conditions
X_{1} = memZono(X_0,xDims(1));
X_all = X_{1};

% Terminal Set
X_F = memZono(X_F,xDims(N));

switch 'fullStack' %'withOverload' 'transform&memoryIntersection' 'fullStack'
    case 'withOverload'
        % Time-evolution
        for k = 1:N-1
            % Current Input
            U_{k} = memZono(U_nom,uDims(k));

            % Step Update
            X_{k+1} = X_{k}.map(A,xDims(k),xDims(k+1)) + U_{k}.map(B,uDims(k),xDims(k+1));

            % Save Data
            X_all = [
                X_all; 
                U_{k}; 
                X_{k+1}]; %<---- vertcat() = cartprod()
        end
        X_inter = X_all.and(X_F,'terminal_cons'); %<-- memoryIntersection does the intersection

    % case 'transform&memoryIntersection'
    %     lbl_ = @(x,k) sprintf('%s_%d',x,k);
    %     % Time-evolution
    %     for k = 1:N-1
    %         % Current Input
    %         % U_{k} = memZono(U_nom,lbl_('u',k));
    %         U_{k} = memZono(U_nom,sprintf('u_%d',k));
    % 
    %         % Step Update
    %         newDims = memZono.genKeys(lbl_('x',k+1),1:n);
    %         % newDims = {sprintf('x_%d_1',k+1),sprintf('x_%d_2',k+1)};
    %         X_{k+1} = X_{k}.transform(...
    %             U_{k}.transform([],B,U_{k}.dimKeys,newDims),A,...
    %             X_{k}.dimKeys,newDims); %<== transform has affine A x + B
    % 
    %         % Save Data
    %         X_all = X_all.memoryIntersection(U_{k});
    %         X_all = X_all.memoryIntersection(X_{k+1});
    %     end
    %     X_inter = X_all.memoryIntersection(X_F,'terminal_cons'); % <--- intersect common dimensions

    case 'fullStack'
        % Time-Evolution
        for k = 1:N-1
            % Current Input
            % X_all = cartProd(X_all,memZono(U_nom,uDims(k)));
            U_k = memZono(U_nom,uDims(k));
            X_all = cartProd(X_all, U_k);
            % X_all = [X_all; U_k];

            % map(), plus(), cartProd()
            % X_new = X_all.map(A,xDims(k),xDims(k+1)) ...
            %         + X_all.map(B,uDims(k),xDims(k+1));
            % X_all = cartProd(X_all,X_new);
            X_all = cartProd(X_all,...
                plus(map(X_all,A,xDims(k),xDims(k+1)),...
                    map(X_all,B,uDims(k),xDims(k+1))));
            % X_all = [X_all;
            %     plus(map(X_all,A,xDims(k),xDims(k+1)),...
            %         map(X_all,B,uDims(k),xDims(k)));
            %     ];
        end

        X_inter = and(X_all,X_F,'terminal_cons');

        X_ = arrayfun(@(k) X_all(xDims(k)),1:N,UniformOutput=false);
        U_ = arrayfun(@(k) X_all(uDims(k)),1:N-1,UniformOutput=false);

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