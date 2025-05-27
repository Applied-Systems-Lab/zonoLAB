%% memZono Example: SLAM
clear; clc;
%% Plot settings 
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultlegendinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')
rng(1); %<== reset random number generator

%  Simulation Settings
N = 5;         % number of time steps (starting at 1)
n_L = 20;     % number of Lnom (must be even)
% theta = pi/6;   % rad between each measurement
% theta = 0.4*pi;

% System Defintion
% A = [cos(theta),sin(theta);-sin(theta),cos(theta)]; % Rotation Matrix
A = eye(2); xshift = [2;0];
B = eye(2); % system noise input
C = eye(2); % Measurement
n = size(A,1); m = size(B,2); p = size(C,1);

% Noise
eta_0 = 1e-1;
V0 = zono(eta_0.*[eye(n),ones(n,1),[-1;1]],zeros(n,1));
nu_0 = 1e-1.*[5,1];
% 5e0*[0.075, 0.05];
H0 = zono(diag(nu_0),zeros(p,1));

% Initial Conditons
x0 = [0;0]; % initial state
X0 = x0 + 1e-1*V0; % initial state set

% Landmarks
L{N,n_L} = [];
r_detection = 5;
% for i = 1:n_L
%     theta = 2*pi*rand();
%     Lnom{i} = (2+2*rand())*[cos(theta);sin(theta)];
% end
lb = [-1;-5]; ub = [10;5];
% Lmeas_0 = H_0;%zono(diag(eta_0),zeros(p,1));
% for i = 1:n_L
%     Lnom{i} = ((ub-lb).*rand(n,1)+lb); %<== random location
%     % while abs(Lnom{i}(2)) <= ub(2)./4
%     %     % Lnom{i} = Lnom{i}.*4;
%     %     Lnom{i} = ((ub-lb).*rand(n,1)+lb); %<== random location
%     % end
% end

Lx = linspace(lb(1),ub(1),n_L/2);
% Ly = linspace(lb(2),ub(2),n_L);
Lbase = [[Lx,Lx];[lb(2).*ones(1,n_L/2)/2, ub(2).*ones(1,n_L/2)/2]];
for i = 1:n_L
    Lnom{i} = Lbase(:,i) + (2*rand(n,1)-1).*[0;1]; %<= normal about the places...
end


% Labeling
xLbl = @(k) sprintf('x_%dk',k);
xDim = @(k) memZono.genKeys(xLbl(k),1:n);
vLbl = @(k) sprintf('v_%dk',k);
vDim = @(k) memZono.genKeys(vLbl(k),1:n);
vKey = @(k) sprintf('v_%dk',k);
lDim = @(i) sprintf('L_%di',i);
lKey = @(k,i) sprintf('L_%dk_%di',k,i);
lCon = @(k,i) sprintf('conL_%dk_%di',k,i);

%% Simulation
X{1} = memZono(X0,xDim(1),xLbl(1)); 
Xall = X{1}; 
x(:,1) = X0.c;
% Time Evolution
for k=1:N %(evolve past N to get measurements)%
    % Landmarks
    for i = 1:n_L
        r = Lnom{i} - x(:,k);
        if norm(r) <= r_detection
            theta{k,i} = atan2(r(2),r(1));
            RH{k,i} = rotate_zonotope(H0,theta{k,i}); % Rotate measuremet
            r_m{k,i} = r + random_sample_zonotope(RH{k,i}); % <==noisy measurment realization
            % Predict based on measurements
            L{k,i} = relabelDims(Xall,xDim(k),lDim(i)) ... %<= relative to state position
                + memZono(RH{k,i}+r_m{k,i},lDim(i),lKey(k,i)); %<= measured distance (r_m) + rotated error measurement (RH)
            % Include landmark
            Xall = Xall.and(L{k,i}, lCon(k,i));
        end
    end 
    Xall_{k} = Xall; %<== save data snapshot

    % Time-Update
    V{k} = memZono(V0,vDim(k),vLbl(k)); % System Noise at k
    % X{k+1} = X{k}.map(A,xDim(k),xDim(k+1)) ... %<= rotate/new position
    X{k+1} = X{k}.map(A,xDim(k),xDim(k+1)) ... %<= new position
            + memZono(xshift,xDim(k+1),'none') ...
            + V{k}.map(B,vDim(k),xDim(k+1)); % <= System Noise Step-update
    Xall = [Xall; X{k+1}];
    x(:,k+1) = A*x(:,k) + B*random_sample_zonotope(V0) + xshift; %<== actual realization
end


clr = flipud(prism);
% clr = jet(N);
% clr = hsv(N);

%% Ploting
fig = figure();
t = tiledlayout('flow');
for k = 1:N
    ax{k} = nexttile; hold on;
    % title(sprintf('$k = %i$',k-1),'interpreter','Latex');
    xlabel('x');
    ylabel('y');
    legend();
    yline(0,'--','HandleVisibility','off');
    grid on;
    axis equal;


    % Nominal
    for i = 1:k-1
       copyobj([p_Xnom_{i}],ax{k});%<== don't want plotted after any others
    end
    [p_Xnom_{k},p_xnom_{k}] = plotZonoAndCentroid(X{k},xDim(k),'k',0.3,'.');
    p_xnom_{k}.Visible = 'off'; %<= don't show centers

    % Copy past estimates
    for i = 1:k-1
        copyobj([p_X_{i,:}],ax{k});
    end
    
    % Plot State Zonotopes
    for j = 1:k
        [p_X_{k,j},p_x_{k,j}] = plotZonoAndCentroid(Xall_{k},xDim(j),clr(k,:),0.3,'+');
        p_x_{k,j}.Visible = 'off'; %<== don't show estimated value
        if j == k
            p_X_{k,j}.FaceAlpha = 0.9;
            p_X_{k,j}.LineWidth = 0.5;
            p_X_{k,j}.DisplayName = sprintf('$k=%d$',k);
            p_X_{k,j}.HandleVisibility = 'on';
        end
    end
    % plot actual
    p_xactual_{k} = plot(x(1,1:k),x(2,1:k),'kx',DisplayName = 'Actual');
    drawnow;
    
    % Plot landmark measuremetns
    % past measurement and estimates
    for j = 1:k-1
        copyobj([p_Lmeas_{j,:}],ax{k}); % measurments at past timestep
        copyobj([p_Lest_{j,:}],ax{k}); % estimates at past timestep
    end
    for i = 1:n_L % Loop Landmarks
        k_first = find(~cellfun(@isempty,L(1:k,i)),1,"first"); 
        k_last = find(~cellfun(@isempty,L(1:k,i)),1,"last");
        if ~isempty(k_last)
            % Plot landmark measurements
            if k_last == k
                [p_Lmeas_{k,i},P_lmeas_{k,i}] = plotZonoAndCentroid(L{k,i},lDim(i),clr(k,:),0,'.');
                p_lmeas_{k,i}.Visible = 'off'; %<== don't show centroid value
            end
        end
    end
    drawnow;

    % Landmark Estimates
    for i = 1:n_L % Plot for k
        if Xall_{k}.Z(lDim(i)).n > 0
            %% plot landmark estimate
            [p_Lest_{k,i},p_lest_{k,i}] = plotZonoAndCentroid(Xall_{k},lDim(i),clr(k,:),0.5,'+');
            p_Lest_{k,i}.EdgeColor='k';
            p_lest_{k,i}.Visible = 'off'; %<== don't show estimated value    
            % Actual Location
            p_lnom_{k,i} = plot(Lnom{i}(1),Lnom{i}(2),'kx',HandleVisibility='off'); % plot Lnom accurate location
        end
    end
    % Plot past landmark estimates
    for j = 1:k; copyobj([p_lest_{j,:}],ax{k}); end
    % Put past postions on top
    for i = 1:k; copyobj([p_x_{i}],ax{k}); end
    drawnow;
end
fig2 = figure; t2 = tiledlayout('flow','TileSpacing','compact','Padding','compact');
copyobj(ax{end},t2);
legend();
title([])
xlabel('x')
ylabel('y')
xlim([-2 11])
ylim([-5 5])
drawnow;

%%% Save Fig
saveas(fig2,'ex_SLAM.png')
% % print to eps file:
% print(fig2,'ex_SLAM','-depsc','-vector')
% exportgraphics(fig2,'ex_SLAM_exportGraphics.eps', ContentType='vector')


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