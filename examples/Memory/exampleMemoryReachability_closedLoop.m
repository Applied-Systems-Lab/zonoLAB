% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example: 
%   Closed-loop Reachability with uncertain values being within sets
% 
% System:
%     dx = A x + B u + \nu, \nu \in V
%     y = C x + \eta, \eta \in H
%     xhat = A xhat + B u + L(y - C xhat)
%     u = K xhat
% 
%     x(t=0) = x_0, x_0 \in X_0
% 
% Reachability:
%     Each set propogates foward... (insert more details)
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Simulation Settings
tf = 1;
dt = 0.02;
tspan = 0:dt:tf;
N = length(tspan);

%% System Definition
m = 1; %Kg
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

% Observer/Controller Definition
L = place(A',C',[0.8;0.85]); 
K = -place(A,B,[0.8;0.5]);

%% Set definition
% Initial Position
X0 = conZono;
X0.c = [1; 0];
X0.G = zeros(n,1);
Xhat0 = X0;
Xhat0.G = 1e-3*eye(n);

% System Noise
V = conZono;
V.G = 2e-3*[0.25, 0.5; 1, -0.25];
V.c = -0.25e-3*ones(size(V.G,1),1);

% Measurment Noise
mu(1) = -1e-3; sigma(1) = 1e-1;
mu(2) = 1e-2; sigma(2) = 5;
H = dt * conZono(diag(sigma),mu');

% Time varying
[V_{N}, H_{N}] = deal([]);
for k = 1:N
    V_{k} = conZonoM(V,sprintf('V_k=%d',k));
    H_{k} = conZonoM(H,sprintf('H_k=%d',k));
end

% Simulation
[X_{N+1}, Xhat_{N+1}, U_{N}, Y_{N}] = deal([]);
X_{1} = conZonoM(X0,'X_k=0'); 
Xhat_{1} = conZonoM(Xhat0, 'Xhat_k=0');
for k = 1:N
    Y_{k} = C*X_{k} + H_{k};
    U_{k} = K*Xhat_{k};
    
    X_{k+1} = A*X_{k} + B*U_{k} + V_{k};
    Xhat_{k+1} = A*Xhat_{k} + B*U_{k} ...
        + L*(Y_{k} + -C*Xhat_{k});
end

%% Ploting
figure;
hold on
for k = 1:N
    plot(Xhat_{k}.Z)
end
