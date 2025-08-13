%% memZono Example: Observer
clear; clc;
close all;

rng(1);


%% Simulation Settings
N = 6; %< number of timesteps to do reachability (index starts at 1)
dt = 0.25;

%% System Definition
% 
% DT-LTI System %%%%
% x_{k+1} = A x_{k} + [B 0] w_{k}
% y_{k}   = C x_{k} + [0 D] w_{k}
% 
% x_{0} \in X_{k}, w_{k} \in W_{k}
% 
A = [   1, 0.25;
     -0.5, 0.75];
B = [   0.025;
        0.25]; shift = [3;4];
C = eye(2);
% D = 0.2*[1.5,0; 0,7];%eye(2);%[1; 1]; 
D = 0.2*[1.5,0.25; 0,7];%eye(2);%[1; 1]; 
% tmp = zeros(size(B)); B = [B, zeros(size(D))];  D = [tmp, D]; %<= no corelation between w and v
tmp = B; B = [B,zeros(size(D))]; D = [10*tmp,D];

n = size(A,1); m = size(B,2); p = size(C,1);
% 
% Estimator Definitions %%%
% 
% Estimator 2
% Xhat_{k+1} = (A Xhat_{k} \oplus B W_{k}) \cap (y_{k} \oplus D W_{k})
% 
% Estimator 3
% Xhat_{k+1} = (A Xhat_{k} \staroplus B W_{k}) \starcap (y_{k} \staroplus D W_{k})

%% Reachability Setup
% Setup Labels
% Labeling
xLbl = @(k) sprintf('x_%dk',k); xDim = @(k) memZono.genKeys(xLbl(k),1:n);
wLbl = @(k) sprintf('w_%dk',k); wDim = @(k) memZono.genKeys(wLbl(k),1:m);
yLbl = @(k) sprintf('y_%dk',k); yDim = @(k) memZono.genKeys(yLbl(k),1:p);
xbarLbl = @(k) sprintf('xbar_%dk',k); xbarDim = @(k) memZono.genKeys(xbarLbl(k),1:n); %<= offline
xhatLbl = @(k) sprintf('xhat_%dk',k); xhatDim = @(k) memZono.genKeys(xhatLbl(k),1:n); %<= intersection
xstarLbl = @(k) sprintf('xstar_%dk',k); xstarDim = @(k) memZono.genKeys(xstarLbl(k),1:n); %<= trajectory
wbarLbl = @(k) sprintf('wbar_%dk',k); %wbarDim = @(k) memZono.genKeys(wbarLbl(k),1:m);
whatLbl = @(k) sprintf('what_%dk',k); %whatDim = @(k) memZono.genKeys(whatLbl(k),1:m);

% Initial Conditions
X0 = zono(diag([1,2]),zeros(2,1));
X_{1} = memZono(X0,xDim(1));

% Noise Sets
W = zono(diag([2,2,2]),zeros(m,1));% + [2;0;0];%<= 2 center applies as a shift... making it botentially not fully including 0 anymore...
for k = 1:N
    W_{k} = memZono(W,wDim(k),wLbl(k));
    Wbar_{k} = memZono(W,wDim(k),wbarLbl(k)); %<= fictions/unknown version
    % What_{k} = memZono(W,wDim(k),whatLbl(k));
end


%% Reachability Analysis
%  Startup
X_{1} = memZono(X0,xDim(1),xLbl(1));
% Xbar_{1} = memZono(X0,xbarDim(1),xbarLbl(1));
% Xhat_{1} = memZono(X0,xhatDim(1),xhatLbl(1));
Xstar_{1} = memZono(X0,xstarDim(1),xstarLbl(1));

% Actual reachability rollout
x_{1} = X0.c; %<= actual position
for k = 1:N
    w_{k} = random_sample_zonotope(W_{k}.Z(W_{k}.dimKeys)); % noise
    x_{k+1} = A*x_{k} + B*w_{k} + shift; % time-update
    y_{k} = C*x_{k} + D*w_{k}; % measurement
    Y_{k} = memZono(y_{k},yDim(k),yLbl(k)) ...
        + W_{k}.map(-D,wDim(k),yDim(k)); %<= measurement w/ uncertainty
end

% Estimator 1 (Offline/Pure Reachability)
Xbar_{1} = memZono(X0,xDim(1),xbarLbl(1));
for k = 1:N
    Xbar_{k+1} = Xbar_{k}.map(A,xDim(k),xDim(k+1)) ...
        + memZono(shift,xDim(k+1))...
        + Wbar_{k}.map(B,wDim(k),xDim(k+1)); % time-update
end

% % Estimator 2 (Online Intersection)
% Xhat_{1} = memZono(X0,xDim(1),xhatLbl(1));
% for k = 1:N
%     % estimated measurement
%     Yhat_{k} = Xhat_{k}.map(C,xDim(k),yDim(k));
%     % intersection w/ Y_{k}
%     Yinter_{k} = and(Yhat_{k},Y_{k},yDim(k));
%     % augment the constraints
%     Xhat_{k} = projection(cartProd(Xhat_{k},Yinter_{k}),xDim(k)); %<== hopefully good enough to do generalized intersection...
    
%     % prediction update
%     Xhat_{k+1} = Xhat_{k}.map(A,xDim(k),xDim(k+1)) ...
%         + memZono(shift,xDim(k+1))...
%         + What_{k}.map(B,wDim(k),xDim(k+1)); %<= hat means no memory w/ actual
% end

Xhat_{1} = and(X0,y_{1}+-1*D*W,C);
for k = 1:N-1
    % % prediction update
    % Xhat_{k+1} = A*Xhat_{k} + shift + B*W;
    % % intersection w/ Y_{k}
    % Xhat_{k+1} = and(Xhat_{k+1}, y_{k+1} + -1*D*W,C); 
    Xhat_{k+1} = and(A*Xhat_{k} + shift + B*W, y_{k+1} + -1*D*W,C); 
end

for k = 1:N
    Xhat_{k} = memZono(Xhat_{k},xDim(k),xLbl(k));
end


% Estimator 3 (Memory)
Xstar_{1} = memZono(X0,xDim(1),xstarLbl(1));
for k = 1:N
    % estimated measurement
    Ystar_{k} = Xstar_{k}.map(C,xDim(k),yDim(k));
    % intersection w/ Y_{k}
    Ystarinter_{k} = and(Ystar_{k},Y_{k},yDim(k));
    % augment the constraints
    Xstar_{k} = projection(cartProd(Xstar_{k},Ystarinter_{k}),xDim(k)); %<== hopefully good enough to do generalized intersection...
    
    % prediction update
    Xstar_{k+1} = Xstar_{k}.map(A,xDim(k),xDim(k+1)) ...
        + memZono(shift,xDim(k+1))...
        + W_{k}.map(B,wDim(k),xDim(k+1));
end


%% Ploting
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultlegendinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

fig = figure();
t = tiledlayout('flow');
for k = 1:N
    ax{k} = nexttile; hold on;
    % title(sprintf('$k = %i$',k-1),'interpreter','Latex');
    xlabel('x');
    ylabel('y');
    % legend();
    yline(0,'--','HandleVisibility','off');
    grid on;
    axis equal;

    % Measurments
    for i = 1:k-1; copyobj([p_y_{i},p_Y_{i}],ax{k}); end
    [p_Y_{k},p_y_{k}] = plotZonoAndCentroid(Y_{k},yDim(k),'b',0.3,'*');


    % Estimator 1
    for i = 1:k-1; copyobj([p_xbar_{i},p_Xbar_{i}],ax{k}); end
    [p_Xbar_{k},p_xbar_{k}] = plotZonoAndCentroid(Xbar_{k},xDim(k),'y',0.3,'.');
        % p_Xbar_{k}.HandleVisibility = "on";
        % p_Xhat_{k}.DisplayName = sprintf('$\overline{\mathcal{X}}_{%d}$',k);

    % Estimator 2
    for i = 1:k-1; copyobj([p_xhat_{i},p_Xhat_{i}],ax{k}); end
    [p_Xhat_{k},p_xhat_{k}] = plotZonoAndCentroid(Xhat_{k},xDim(k),'r',0.3,'*');
        % p_Xhat_{k}.HandleVisibility = "on";
        % p_Xhat_{k}.DisplayName = sprintf('$\hat{\mathcal{X}}_{%d}$',k);

    % Estimator 3
    for i = 1:k-1; copyobj([p_xstar_{i},p_Xstar_{i}],ax{k}); end
    [p_Xstar_{k},p_xstar_{k}] = plotZonoAndCentroid(Xstar_{k},xDim(k),'g',0.3,'o');
        % p_Xstar_{k}.HandleVisibility = "on";
        % p_Xstar_{k}.DisplayName = sprintf('$\tilde{\mathcal{X}}_{%d}$',k);


    % Draw Actual Locations
    for i = 1:k; plot(x_{i}(1),x_{i}(2),'kx'); end

    drawnow;
end
fig2 = figure; t2 = tiledlayout('flow','TileSpacing','compact','Padding','compact');
copyobj(ax{end},t2);
% legend();
title([])
xlabel('x')
ylabel('y')
axis square
% xlim([-2 11])
% ylim([-5 5])
drawnow;





%% Local Functions
function Z_prime = rotate_zonotope(Z, radians)
    center = Z.c;
    Z.c = Z.c - center;
    R = [cos(radians) -sin(radians); sin(radians) cos(radians)];
    Z_prime = R * Z;
    Z_prime.c = Z_prime.c + center;
end

% random_sample_zonotope()
function s = random_sample_zonotope(z)
    g = 2*rand([z.nG,1])-1;
    s = z.c + z.G*g;
end


function [P,p] = plotZonoAndCentroid(Zm,dim,clr,w,mkr,varargin)
    [v,~] = plot(Zm,dim,clr,w);
    ax = gca; P = ax.Children(1);
    P.HandleVisibility='off';
    P.EdgeAlpha=1;
    if w == 0; P.EdgeColor=clr; else; P.EdgeColor = 'k'; end
    [x,y] = centroid(polyshape(v));
    p = plot(x,y,mkr,...
        "MarkerEdgeColor",clr,...
        "HandleVisibility","off",varargin{:});
end