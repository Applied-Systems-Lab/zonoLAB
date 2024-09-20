%% memZono Example: SLAM
clear; clc;
rng(1); %<== reset random number generator

%% Setup
%  Simulation Settings
N = 4;         % number of time steps (starting at 1)
n_L = 30;     % number of Lnom
% theta = pi/6;   % rad between each measurement
theta = 0.4*pi;

% System Defintion
A = [cos(theta),sin(theta);-sin(theta),cos(theta)]; % Rotation Matrix
B = eye(2); % system noise input
C = eye(2); % Measurement
n = size(A,1); m = size(B,2); p = size(C,1);

% Noise
eta_0 = 0.05;
V0 = zono(eta_0*eye(n),zeros(n,1));
nu_0 = [0.075, 0.05];
H0 = zono(diag(nu_0),zeros(p,1));

% Initial Conditons
x0 = [5;0]; % initial state
X0 = zono(1e-2*eye(n),x0); % initial state set

% Landmarks
L{N,n_L} = [];
r_detection = 4;
for i = 1:n_L
    theta = 2*pi*rand();
    Lnom{i} = (2+2*rand())*[cos(theta);sin(theta)];
end

% Labeling
xDim = @(k) memZono.genKeys(sprintf('x_%dk',k),1:n);
vDim = @(k) memZono.genKeys(sprintf('v_%dk',k),1:n);
vKey = @(k) sprintf('v_%dk',k);
lDim = @(i) sprintf('L_%di',i);
lKey = @(k,i) sprintf('L_%dk_%di',k,i);
lCon = @(k,i) sprintf('conL_%dk_%di',k,i);

%% Simulation
X{1} = memZono(X0,xDim(1)); 
Xall = X{1}; 
x(:,1) = X0.c;
% Time Evolution
for k=1:N %(evolve past N to get measurements)%
    % Landmarks
    for i = 1:n_L
        r = Lnom{i} - x(:,k);
        if norm(r) <= r_detection
            theta = atan2(r(2),r(1));
            RH = rotate_zonotope(H0,theta); % Rotate measuremet
            r_m{k,i} = r + random_sample_zonotope(RH); % <==noisy measurment realization
            % Predict based on measurements
            L{k,i} = relabelDims(X{k},xDim(k),lDim(i)) ... %<= relative to state position
                + memZono(RH+r_m{k,i},lDim(i),lKey(k,i)); %<= measured distance (r_m) + rotated error measurement (RH)
            % Include landmark
            Xall = Xall.and(L{k,i}, lCon(k,i));
        end
    end 
    Xall_{k} = Xall; %<== save data snapshot

    % Time-Update
    V{k} = memZono(V0,vDim(k)); % System Noise at k
    X{k+1} = X{k}.map(A,xDim(k),xDim(k+1)) ... %<= rotate/new position
            + V{k}.map(B,vDim(k),xDim(k+1)); % <= System Noise Step-update
    Xall = [Xall; X{k+1}];
    x(:,k+1) = A*x(:,k) + B*random_sample_zonotope(V0); %<== actual realization
end



%% Ploting
fig = figure();
t = tiledlayout('flow');
for k = 1:N
    ax = nexttile; hold on;
    title(sprintf('$k = %i$',k-1),'interpreter','Latex');
    xlabel('x');
    ylabel('y');
    
    % Loop Landmarks
    for i = 1:n_L
        k_last = find(~cellfun(@isempty,L(1:k,i)),1,"last");
        if ~isempty(k_last)
            if k_last == k
                plot(r_m{k,i}(1)+x(1,k),r_m{k,i}(2)+x(2,k),'b.','MarkerSize',8);
                plot([x(1,k),r_m{k,i}(1)+x(1,k)],[x(2,k),r_m{k,i}(2)+x(2,k)],'b','LineWidth',0.5);
                drawnow;
            end
            plot(Xall_{k},lDim(i),'r',0.6);    % plot landmark zonotope
            plot(Lnom{i}(1),Lnom{i}(2),'bx'); % plot Lnom accurate location
            drawnow;
        end
    end

    % Plot State Zonotopes
    for i = 1:k
        plot(x(1,i),x(2,i),'k','MarkerSize',12);
        plot(X{i},xDim(i),'k',0.2);
        plot(Xall_{k},xDim(i),'g',0.3);
    end

end





%%% Save Fig:
saveas(fig,'ex_SLAM.png')






%% Plot Function
% function plotSLAM(Xall,k,i)












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