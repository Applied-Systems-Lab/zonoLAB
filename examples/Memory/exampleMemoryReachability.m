% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example: 
%   Reachability with uncertain values being within sets
% 
% System:
%     dx = A x + B u + \nu, \nu \in V
%     y = C x + \eta, \eta \in H
% 
% Reachability:
%     Each set propogates foward... (insert more details)
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear

%% Simulation Settings
tf = 3;
dt = 0.02;
tspan = 0:dt:tf;
N = length(tspan);

%% System Definition
m = 0.1; %Kg
k = 1;   %N/m
b = 0.25; %N/m^2
A_ct = [0, 1; -k/m, -b/m];
B_ct = [0; -1/m];
C = eye(2); D = 0;
sys_ct = ss(A_ct,B_ct,C,D);

% DT conversion
sys = c2d(sys_ct,dt);
A = sys.A; B = sys.B;

% Dims
n = size(A,1); 
m = size(B,2); 
p = size(C,1);

%% Set definition
% Initial Position
X0 = conZono;
X0.c = [0.5; 0];
X0.G = 1e-2*[1; 0];%1e-2*eye(n);
% Xhat0 = X0;
% Xhat0.G = 1e-3*eye(n);

% % System Noise
% V = conZono;
% V.G = 2e-3*[0.25, 0.5; 1, -0.25];
% V.c = -0.25e-3*ones(size(V.G,1),1);
% 
% % Measurment Noise
% mu(1) = -1e-3; sigma(1) = 1e-1;
% mu(2) = 1e-2; sigma(2) = 5;
% H = dt * conZono(diag(sigma),mu');

% Input Range
U = conZono;
U.c = zeros(m,1);
U.G = 1e-2*eye(m);

% Time varying
% [V_{N}, H_{N}, U_{N}] = deal([]);
U_{N} = [];
for k = 1:N
    % V_{k} = conZonoM(V,sprintf('V_k=%d',k));
    % H_{k} = conZonoM(H,sprintf('H_k=%d',k));
    U_{k} = conZonoM(U,sprintf('U_k=%d',k));
    % U_{k}.dimKeys = sprintf('u_k=%d',k);
end

% Simulation
% [X_{N+1}, Y_{N}] = deal([]);
X_{N+1} = [];
X_{1} = conZonoM(X0,'X0');
for k = 1:N
    % Y_{k} = C*X_{k} + H_{k};
    
    % X_{k+1} = A*X_{k} + B*U_{k} + V_{k};
    % X_{k+1} = A*X_{k} + V_{k};
    X_{k+1} = A*X_{k} + B*U_{k};

    % label specific rows
    X_{k}.dimKeys = sprintf('x_k=%d',k);
    U_{k}.dimKeys = sprintf('u_k=%d',k);
end
X_{k+1}.dimKeys = sprintf('x_k=%d',k+1);

% Full Info
Xall = vertcat(X_{:}); %<-- cart prod
Uall = vertcat(U_{:}); %<-- cart prod
data = [Xall; Uall];

%% Final intersection
Xf = X0; Xf.c = zeros(n,1);
Xf = conZonoM(Xf,'Xf');
Xf.dimKeys = sprintf('x_k=%d',N);
% for i = 1:n; dims{i} = sprintf('x_k=%d_%d',N,i); end
data_XfInt = labeledIntersection(data,Xf,Xf.dimKeys,'Xf_int');



%% Plotting
figure;
hold on
for k = 1:5:N
    plot(X_{k}.Z,'k')
    % add dims based plot?
    Xk = data_XfInt({...
        sprintf('x_k=%d_%d',k,1),...
        sprintf('x_k=%d_%d',k,2)});
    plot(Xk.Z,'b')
end
