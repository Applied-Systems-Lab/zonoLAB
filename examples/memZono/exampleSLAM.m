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
rng(5); %<== reset random number generator

% %% Setup
% N = 5; %<= number of timesteps
% n = 2; %<= 2D vs 3D
% % m = n; % <== input (none)
% p = n; % <== measurement

% % Labeling
% lblDim = @(lbl,n) memZono.genKeys(lbl,1:n);
% % states (x)
% xLbl = @(k) sprintf('x_%dk',k);
% xDim = @(k) lblDim(xLbl(k),n);
% % system noise (nu)
% nuLbl = @(k) sprintf('nu_%dk',k);
% nuDim = @(k) lblDim(nuLbl(k),n);
% % measurement noise (eta) <= position...
% etaLbl = @(k,i) sprintf('eta_%dk_%dl',k);
% etaDim = @(k) lblDim(sprintf('eta_%d',k),n);
% % landmarks (l)
% lLbl = @(k,i) sprintf('l_%dk_%di',k,i);
% lDim = @(i) lblDim(sprintf('l_%di',i),n);

% %% Noise
% % system noise (nu)
% nu_0 = 0*[2;0.5];
% V_0 = zono(nu_0.*[eye(n),ones(n,1),[1;-1]],zeros(n,1)); 
% for k = 1:N-1
%     nu_{k}=nu_0.*(2*rand(n,1)-1);
%     V_{k} = memZono(nu_{k}+V_0,nuDim(k),nuLbl(k)); 
% end
% % measurement noise (eta) <= landmark position...
% eta_0 = 1e-1*ones(n,1);
% H_0 = zono(eta_0.*[eye(n),ones(n,1),[1;-1]],zeros(n,1));
% % for k = 1:N
% %     eta_{k}=eta_0.*(2*rand(n,1)-1);
% %     H_{k} = memZono(eta_{k}+H_0,etaDim(k),etaLbl(k));
% % end

% %% Landmark Setup
% n_L = 20;
% r_detection = 5;
% lb = [-7;-5]; ub = [25;5];
% % Lmeas_0 = H_0;%zono(diag(eta_0),zeros(p,1));
% for i = 1:n_L
%     lnom_{i} = ((ub-lb).*rand(n,1)+lb); %<== random location
%     % for k = 1:N
%     %     lmeas_{k,i} = lnom_{i} + eta_0.*rand(n,1); %<== random error measurement
%     %     H_{k,i} = memZono(H_0,etaDim(k),etaLbl(k,i));
%     %     Lmeas_{k,i} = memZono(lmeas_{k,i},lDim(i),lLbl(k,i)) ...
%     %             + relabelDims(H_{k},etaDim(k),lDim(i));
%     % end
% end

% % initial conditions
% % xstep
% xstep_0 = [5;0];
% for k=1:N-1; xstep_{k} = xstep_0; end
% %xnom
% xnom_0 = zeros(n,1); 
% xnom_{1} = xnom_0 + xstep_0;
% for k = 1:N-1; xnom_{k+1} = xnom_{k} + xstep_{k} + sum([xstep_{1:k}],2); end
% %Xnom
% Xnom_0 = memZono(zono(xnom_0) + V_0, xDim(0),xLbl(0));
% Xnom_{1} = memZono(xnom_{1} + V_0, xDim(1),xLbl(1));
% % xactual
% xactual_{1} = xnom_0;% + xstep_0 + nu_0.*(2*rand(n,1)-1);
% %Zm
% Zm = Xnom_0; %<== only state at 0

% % Object Position
% for k=1:N-1
%     Xstep_{k} = memZono(xstep_{k},xDim(k+1));%<= k+1 as it goes from k -> k+1% state update
%     xactual_{k+1} = xactual_{k} + xstep_{k};% + nu_{k};
%     % Xnom_{k+1} = Xnom_{k}.relabelDims(xDim(k),xDim(k+1)) ...
%     %         + Xstep_{k} + V_{k}.relabelDims(nuDim(k),xDim(k+1));
% end
% % for k = 1:N
% %     Xactual_{k} = memZono(xactual_{k},xDim(k),[xLbl(k),'_actual']);
% %     % Xnom_{k} = Xactual_{k} + relabelDims(-1*memZono.sum(H_(k,:)),etaDim(k),xDim(k));
% % end

% % Time Evolution for SLAM
% L_{N,n_L} = [];
% for k = 1:N
%     % % Add nominal estimate
%     % Zm = cartProd(Zm,Xnom_{k});
%     % Measure landmarks
%     for i = 1:n_L
%         if norm(lnom_{i}-xactual_{k}) < r_detection;
%             eta_k = eta_0.*rand(n,1);
%             lmeas_{k,i} = lnom_{i} + eta_k; %<== random error measurement
%             H_{k,i} = memZono(H_0,etaDim(k),etaLbl(k,i));
%             L_{k,i} = memZono(lmeas_{k,i},lDim(i),lLbl(k,i)) ...
%                     + relabelDims(H_{k,i},etaDim(k),lDim(i));
%             % L_{k,i} = Lmeas_{k,i};%memZono(Lmeas_{k,i},lDim(i),lLbl(k,i)); %<= L_{k,i} is estimaated set
%             % \starcap operation
%             Zm = and(Zm,L_{k,i},lLbl(k,i));
%             % use measurement to update
%             Xest_k = memZono(xactual_{k} - eta_k,xDim(k),xLbl(k)) ...
%                     + relabelDims(H_{k,i},etaDim(k),xDim(i));
%             Zm = and(Zm,Xest_k,lLbl(k,n_L+i));
%         end
%     end
%     Zm_{k} = Zm;
% end




% plot_k = 1:N;

% figure;
% t = tiledlayout('flow');
% for k = plot_k
%     % create new plot
%     ax{k} = nexttile; hold on;
%     xlim([lb(1) ub(1)]);
%     ylim([lb(2) ub(2)]);
%     if k>1
%         % copyobj(ax{k-1}.Children,ax{k}); 
%         p_Xold = copyobj([p_Xest{:}],ax{k});
%         p_Lold = copyobj([p_Lest{:}],ax{k});

%         set([p_Xold,p_Lold],'FaceColor','k');
%     end
%     % Plot Position
%     % Xnom_{k}.plot(xDim(k),'k',0.2); p_Xnom{k} = ax{k}.Children(1);
%     Zm_{k}(xDim(k)).plot(xDim(k),'g',0.2); p_Xest{k} = ax{k}.Children(1);


%     % Plot landmarks
%     for i = 1:n_L
%         if ~isempty(L_{k,i})
%             L_{k,i}.plot(lDim(i),'k',0.2);
%             Zm_{k}.plot(lDim(i),'g',0.1); p_Lest{k} = ax{k}.Children(1);
%             l{k,i} = scatter(lmeas_{k,i}(1),lmeas_{k,i}(2),'xr');
%             lnom{i} = scatter(lnom_{i}(1),lnom_{i}(2),'.k');
%         end
%     end
% end

% fig = figure;
% copyobj(ax{end},gcf)



% return;














%  Simulation Settings
N = 4;         % number of time steps (starting at 1)
n_L = 20;     % number of Lnom (must be even)
% theta = pi/6;   % rad between each measurement
% theta = 0.4*pi;

% System Defintion
% A = [cos(theta),sin(theta);-sin(theta),cos(theta)]; % Rotation Matrix
A = eye(2); xshift = [1;0];
B = eye(2); % system noise input
C = eye(2); % Measurement
n = size(A,1); m = size(B,2); p = size(C,1);

% Noise
eta_0 = 1e-1;
V0 = zono(eta_0*eye(n),zeros(n,1));
nu_0 = 1e-1.*[5,1];
% 5e0*[0.075, 0.05];
H0 = zono(diag(nu_0),zeros(p,1));

% Initial Conditons
x0 = [0;0]; % initial state
X0 = x0 + V0; % initial state set

% Landmarks
L{N,n_L} = [];
r_detection = 3;
% for i = 1:n_L
%     theta = 2*pi*rand();
%     Lnom{i} = (2+2*rand())*[cos(theta);sin(theta)];
% end
lb = [-1;-2]; ub = [5;2];
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
    V{k} = memZono(V0,vDim(k)); % System Noise at k
    % X{k+1} = X{k}.map(A,xDim(k),xDim(k+1)) ... %<= rotate/new position
    X{k+1} = X{k}.map(A,xDim(k),xDim(k+1)) ... %<= rotate/new position
            + memZono(xshift,xDim(k+1),'none') ...
            + V{k}.map(B,vDim(k),xDim(k+1)); % <= System Noise Step-update
    Xall = [Xall; X{k+1}];
    x(:,k+1) = A*x(:,k) + B*random_sample_zonotope(V0) +xshift; %<== actual realization
end


clr = flipud(prism);
% clr = jet(N);
% clr = hsv(N);

%% Ploting
fig = figure();
t = tiledlayout('flow');
for k = 1:N
    ax{k} = nexttile; hold on;
    for i = 1:k-1
       copyobj([p_Xnom_{i}],ax{k});%<== don't want plotted after any others
    end
    for i = 1:k-1
        copyobj([p_X_{i,:}],ax{k});
        copyobj([p_Lest_{i,:}],ax{k});
        copyobj([p_Lmeas_{i,:}],ax{k});
    end
    % if k>1
    %     p_Xnom_old = copyobj(flipud([p_Xnom_{:}]),ax{k});
    %     p_X_old = copyobj(flipud([p_X_{:}]),ax{k});
    %     p_x_old = copyobj(flipud([p_x_{:}]),ax{k});
    %     % % try
    %         p_L_old = copyobj(flipud([p_L_{:}]),ax{k});
    %         p_l_old = copyobj(flipud([p_l_{:}]),ax{k});
    %     %     % set(p_L_old, 'FaceColor','k');
    %     % % catch;end
    %     % % set([p_Xnom_old;p_x_old;p_L_old;p_l_old],'HandleVisibility','off');
    %     % copyobj(ax{k-1}.Children,ax{k});
    % end
    
    % axes setup
    % title(sprintf('$k = %i$',k-1),'interpreter','Latex');
    xlabel('x');
    ylabel('y');
    legend();
    yline(0,'--','HandleVisibility','off');
    grid on;
    
    % Loop Landmarks
    for i = 1:n_L
        k_first = find(~cellfun(@isempty,L(1:k,i)),1,"first"); 
        k_last = find(~cellfun(@isempty,L(1:k,i)),1,"last");
        if ~isempty(k_last)
            % plot(Xall_{k},lDim(i),'r',0.6);    % plot landmark zonotope

            %% Plot landmark measurements
            if k_last == k
                [p_Lmeas_{k,i},P_lmeas_{k,i}] = plotZonoAndCentroid(L{k,i},lDim(i),clr(k,:),0,'.');
                % [v,f] = plot(L{k,i},lDim(i),clr(k,:),0);
                % p_Lmeas_{k,i} = ax{k}.Children(1);
                % p_Lmeas_{k,i}.HandleVisibility = 'off';
                % p_Lmeas_{k,i}.EdgeColor=clr(k,:);
            end
            %% plot landmark estimate
            % if k_first == k; w = 1; elseif k_last == k; w = 0.7; else; w = 0.5; end
            [p_Lest_{k,i},p_lest_{k,i}] = plotZonoAndCentroid(Xall_{k},lDim(i),clr(k,:),0.5,'+');
            p_Lest_{k,i}.EdgeColor='k';
            if k_first == k
                p_Lest_{k,i}.FaceAlpha = 1;
            elseif k_last == k
                p_Lest_{k,i}.FaceAlpha = 0.7;
            end
            % plot(Xall_{k},lDim(i),clr(k,:),w);    % plot landmark zonotope
            % p_Lest_{k,i} = ax{k}.Children(1);
            % p_Lest_{k,i}.HandleVisibility='off';
            % p_lest_{k,i} = plot(Xall_{k}.Z(lDim(i)).c(1),Xall_{k}.Z(lDim(i)).c(2),'+',MarkerEdgeColor=clr(k,:),MarkerSize=12,HandleVisibility='off'); % plot L center
            % plot landmark measurements
            % plot actual
            p_lnom_{k,i} = plot(Lnom{i}(1),Lnom{i}(2),'kx',HandleVisibility='off'); % plot Lnom accurate location
            drawnow;
        end

        % if ~isempty(k_last)
        %     if k_last == k
        %         plot(r_m{k,i}(1)+x(1,k),r_m{k,i}(2)+x(2,k),'b.','MarkerSize',8);
        %         % plot([x(1,k),r_m{k,i}(1)+x(1,k)],[x(2,k),r_m{k,i}(2)+x(2,k)],'b','LineWidth',0.5);
        %         drawnow;
        %     end
        %     plot(Xall_{k},lDim(i),'r',0.6);    % plot landmark zonotope
        %     p_L_{k,i} = ax{k}.Children(1);
        %     plot(Lnom{i}(1),Lnom{i}(2),'bx'); % plot Lnom accurate location
        %     drawnow;
        % end
    end

    % Plot State Zonotopes
    % Nominal (plots only new nominal)
    [p_Xnom_{k},p_xnom_{k}] = plotZonoAndCentroid(X{k},xDim(k),'k',0.3,'x');
    % plot(X{k},xDim(k),'k',0.3);
    % p_Xnom_{k} = ax{k}.Children(1);
    % p_Xnom_{k}.HandleVisibility='off';

    for j = 1:k
        [p_X_{k,j},p_x_{k,j}] = plotZonoAndCentroid(Xall_{k},xDim(j),clr(k,:),0.3,'+');
        % p_x_{k} = plot(x(1,i),x(2,i),'kx','MarkerSize',12,HandleVisibility='off');
        % plot(Xall_{k},xDim(i),clr(k,:),0.3);
        % p_X_{k,i} = ax{k}.Children(1);
        % p_X_{k,i}.HandleVisibility = 'off';
        p_X_{k,j}.EdgeColor = 'k';
        if j == k
            p_X_{k,j}.FaceAlpha = 0.9;
            p_X_{k,j}.LineWidth = 0.5;
            p_X_{k,j}.DisplayName = sprintf('$k=%d$',k);
            p_X_{k,j}.HandleVisibility = 'on';
        end
    end
    % i = k;
    % [p_X_{k},p_x_{k,i}] = plotZonoAndCentroid(Xall_{k},xDim(i),clr(k,:),1,'.');
    % plot(Xall_{k},xDim(i),clr(k,:),1);
    % p_X_{k,i} = ax{k}.Children(1);
    % p_X_{k,k}.DisplayName = sprintf('$k=%d$',k);
    % p_x_{k} = plot(x(1,i),x(2,i),'kx','MarkerSize',12,HandleVisibility='off');

    % Put past postion markers on top
    for i = 1:k-1
        copyobj([p_x_{i}],ax{k});       
        % copyobj([p_lest_{i,:}],ax{k});
    end
end
fig2 = figure; t2 = tiledlayout('flow');
copyobj(ax{end},t2);
legend();
title([])
xlabel("dim = `x'")
xlabel("dim = `y'")


%%% Save Fig:
saveas(fig2,'ex_SLAM.png')






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



function [P,p] = plotZonoAndCentroid(Zm,dim,clr,w,mkr,varargin)
    [v,~] = plot(Zm,dim,clr,w);
    ax = gca; P = ax.Children(1);
    P.HandleVisibility='off';
    if w >= 0.5, P.EdgeColor='k'; else, P.EdgeColor=clr; end;
    P.EdgeAlpha=1;
    [x,y] = centroid(polyshape(v));
    p = plot(x,y,mkr,...
        "MarkerEdgeColor",clr,...
        "HandleVisibility","off",varargin{:});
end