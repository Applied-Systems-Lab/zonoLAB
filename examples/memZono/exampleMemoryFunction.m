%% Example Function Creation
clear; 
% close all;

% %% Function Setup
% n = 11;
% x0 = 1.5*pi;
% bd = [-x0,x0];
% f = @(x) sin(x);
% % f = @(x) atan(x);

% % f_SOS = makeSOSfunc(f,n,bd);
% % F_SOS = memZono(f_SOS,{'x','y'},'sos');
% F_SOS = makeSOSfunc(f,n,bd);
% f_SOS = F_SOS.Z({'x','y'});

% %% Account for error
% errSOS = makeSOSerror(f,n,bd);
% errSOS_max = max(abs(errSOS),[],'all');
% f_Err = zono([0;errSOS_max],zeros(2,1));
% F_Err = memZono(f_Err,{'x','y'},'err');

% %% Construct function
% F = f_SOS+f_Err;
% Fm = F_SOS + F_Err;

% if false
% %% plotting
% figure;
% t = tiledlayout("vertical");
% % f(x)
% dx = 0.1;
% x = bd(1):dx:bd(2);
% y = f(x);
% % sos
% x_sos = linspace(bd(1),bd(2),2*n+1);
% y_sos = f(x_sos);

% %% Plot full region
% ax = nexttile; hold on;
% Fm.plot({'x','y'},'b',0.5);
% plot(x_sos,y_sos, 'k--','DisplayName','f_{SOS}(x)');
% plot(x,y,'r--','DisplayName','f(x)');
% legend


% %% X intersections
% % dx_sos = 0.25;%diff(x_sos([1 2]));
% % for i = 1:numel(x_sos)
% %     X_{i} = zono([dx_sos/4],x_sos(i));
% % end
% % Xall = [X_{:}];
% % X = X_{1};
% % for i = [3,5,7,9,11]
% %     X = union(X,X_{i});
% % end
% % % for i=; X = union(X,Xall(i)); end

% X = zono([pi/2],0);
% Xm = memZono(X,'x','X');

% FXm = and(Fm,Xm,'eval');
% ax = nexttile; hold on;
% plot(FXm,{'x','y'},'r',0.5);
% plot(x,y,'r--');

% % Y intersections
% Y = zono(0.25,0.25);
% Ym = memZono(Y,'y','Y');

% FYm = and(Fm,Ym,'eval');
% ax = nexttile; hold on;
% plot(FYm,{'x','y'},'b',1);
% plot(x,y,'r--');

% %
% Ym2 = memZono(union(Y+0.25,-1*(Y+0.25)),'y','Y2');
% FYm2 = and(Fm,Ym2,'eval');
% ax = nexttile; hold on;
% plot(FYm2,{'x','y'},'b',1);
% plot(x,y,'r--');

% end


% if false
% %% Shifting
% figure; 
% t = tiledlayout('vertical');
% ax = nexttile; hold on;

% %% Domain Shift
% % xshift = pi/2;
% % shift_X = makeFshift(xshift,2*bd);
% % shift_Xm = memZono(shift_X,{'in','out'},'xshift');
% % shift_Xm = makeFshift(xshift,2*bd);
% % shift_X = shiftX.Z({'in','out'});
% % Fm_shiftX = and(...
% %     relabelDims(shift_Xm,{'in','out'},{'x_xshift','inter'}),...
% %     relabelDims(Fm,{'x','y'},{'inter','f_xshift'}),...
% % 'shift');
% xshift = pi/2;
% shift_Xm = makeFshift(xshift,2*bd,{'x','inter'},'xshift');
% Fm_shiftX = and(shift_Xm, ...
%     relabelDims(Fm,{'x','y'},{'inter','y'}), 'shift');


% %% Range shift
% % yshift = 0.5;
% % shift_Y = makeFshift(yshift,[-2 2]);
% % shift_Ym = memZono(shift_Y,{'in','out'},'yshift');
% % shift_Ym = makeFshift(yshift,[-2 2]);
% % shift_Y = shift_Ym.Z({'in','out'});
% % Fm_shiftY = and(...
% %     relabelDims(shift_Ym,{'in','out'},{'inter','f_yshift'}),...
% %     relabelDims(Fm,{'x','y'},{'x_yshift','inter'}),...
% % 'shift');

% yshift = 0.5;
% shift_Ym = makeFshift(yshift,[-2 2],{'inter','y'},'yshift');
% Fm_shiftY = shift_Ym.and( ...
%     relabelDims(Fm,{'x','y'},{'x','inter'}),'shift');

% %% Ploting
% % Fall = [Fm; Fm_shiftX; Fm_shiftY];
% plot(Fm,{'x','y'},'k','0.5');
% plot(Fm_shiftX,{'x','y'},'r',0.5);
% plot(Fm_shiftY,{'x','y'},'b',0.5);


% %% adding F and shift F
% ax = nexttile; hold on;
% plot(Fm,{'x','y'},'k','0.5');
% plot(Fm_shiftX,{'x','y'},'r',0.5);
% Fm_plusShift = Fm + Fm_shiftX;
% Fm_plusShift.plot({'x','y'},'g',0.5);

% end

% if false
% %% Symbolic things
% xshift = 10;
% shiftX = makeFshift(xshift,2*bd,{'x','inter'},'xshift');
% Xsym = memZono(sym('x'),'x','xsym');
% shiftXsym = and(Xsym,shiftX,'shift');
% end


if true

%% Function Setup
n = 11;
x0 = 5;
bd = [-x0,x0];

% functions
f1 = @(x) cos(x);
f2 = @(x) sin(x);
% f = @(x) f1(x(1)) + f2(x(2));

Fm1 = makeSOSfunc(f1,n,bd,{'x1','y'},'sin');
Fm2 = makeSOSfunc(f2,n,bd,{'x2','y'},'cos');
Fm = Fm1 + Fm2;


% NN version
% load("examples\Neural_Networks\relu_sin_cos_2_20_10_10_1.mat");
load("examples\Neural_Networks\NEWsincos_20_10_10.mat")
% load("examples\Neural_Networks\NEW2sincos_20_10_10.mat")
X = memZono(zono(diag([x0,x0]),zeros(2,1)),{'x1','x2'});
NN = reluNN(X,Ws,bs,a);
NN = NN.relabelDims(NN.dimKeys,{'x1','x2','y'}); %<== assumes defaults are in order

% Exaluation Version
xstep = 0.1;
x = -x0:xstep:x0;
[X1,X2] = meshgrid(x,x);
Y = f1(X1) + f2(X2);

% error
Y_Fm = X.and(Fm,'Xin').projection('y'); %<= map X through Fm and project to y
Y_NN = NN('y'); %<= projection of NN along Y (already maps X)
Y_E = Y_NN + -1*Y_Fm; %<== calc difference between them (w/ same input space factors)
E = cartProd(X,Y_E); %<= combine together again...
% end


% if true
% Plot
dims = {'x1','x2','y'};

figure;
t = tiledlayout("flow");
% Actual
ax(1) = nexttile; hold on;
surf(X1,X2,Y, 'EdgeColor', 'none');
title('f(x) = cos(x_1) + sin(x_2)');

% Func Creation
ax(end+1) = nexttile; hold on;
plot(Fm,dims,'r',1);
title('memZono Direct Construction')

% Func error
ax(end+1) = nexttile; hold on;
plot(E,dims,'r');
title('memZono Direct Error')


% NN version
ax(end+1) = nexttile; hold on;
plot(NN.Z(dims),'r',1);
title('trained reluNN construction')

% % NN error
% Y_NN = arrayfun(@(x1,x2) FmEval(NN,[x1;x2]),X1,X2);
% ax(end+1) = nexttile; hold on;
% surf(X1,X2,Y_NN - Y);
% title('NN error')


%% common settings
view(ax,3);
xlabel(ax,dims{1});
ylabel(ax,dims{2});
zlabel(ax,dims{3});
zlim(ax,[-2 2]);

end






if false
    % sin(x) - cos(x)
    n=11; x0 = 5; bd = [-x0,x0];
    F_sin = makeSOSfunc(@sin,n,bd,{'x','y'},'sin');
    F_cos = makeSOSfunc(@cos,n,bd,{'x','y'},'cos');
    Xin = memZono(zono(x0,0),'x');
    Y_sin = Xin.and(F_sin,'xsin').projection('y');
    Y_cos = Xin.and(F_cos,'xcos').projection('y');
    Y_sin_minus_cos = Y_sin + -1*Y_cos;
    XY_cos_sin = [Xin; Y_sin_minus_cos]
    
    Y_vol = memZono(zono(0.1,0),'y');
    
    figure; hold on;
    plot(F_sin+Y_vol,{'x','y'},'r');
    plot(F_cos+Y_vol,{'x','y'},'b');
    plot(XY_cos_sin+Y_vol,{'x','y'},'k');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeSOSfunc
function F_SOS = makeSOSfunc(f,n,bd,dimKeys,lbl)
    arguments
        f
        n 
        bd
        dimKeys = {'x','y'};
        lbl = 'f_sos'
    end
    th = linspace(bd(1), bd(2), 2*n+1);
    xi = th;
    yi = f(th);
    
    V = [xi; yi];
    nc = 2*n+1;
    nb = nc-1;
    naux = nc;
    c = zeros(2,1);
    Gc = V;
    
    M = zeros(nc,nb);
    for i = 1:nb
        M(i,i) = 1;
        M(i+1,i) = 1;
    end
    
    c = c+0.5*Gc*ones(nc,1);
    Gc = 0.5*Gc;
    Gc = [Gc,zeros(2,naux)];
    Gb = zeros(2,nb);
    
    Ac = [0.5*ones(1,nc) zeros(1,naux)]; %(1)
    Ab = zeros(1,nb);
    b = [1-0.5*nc];
    
    Ac = [Ac; 0.5*eye(nc) 0.5*eye(nc)]; %(2)
    Ab = [Ab; -0.5*M];
    b = [b; 0.5*M*ones(nb,1)-ones(nc,1)];
    
    Ac = [Ac; zeros(1,nc+naux)]; %(3)
    Ab = [Ab; 0.5*ones(1,nb)];
    b = [b; 1-0.5*nb];
    
    F_SOS = memZono(hybZono(Gc,Gb,c,Ac,Ab,b),dimKeys,lbl);
end

% makeSOSfunc
function errSOS = makeSOSerror(f,n,bd)
    x = linspace(bd(1),bd(2),n);
    errSOS = zeros(n-1,2);
    for i = 1:n-1
        x0 = x(i);
        x1 = x(i+1);
        x2 = (x0+x1)/2;
        
        y0 = f(x0);
        y1 = f(x1);                
        y2 = f(x2);

        ysos2 = (y0+y1)/2;
        errSOS(i) = y2-ysos2;
    end
end

% x-shift
function F_shift = makeFshift(shift,bnd,dimKeys,lbl)
    arguments
        shift
        bnd = [-10,10]
        dimKeys = {'in','out'};
        lbl = 'fshift';
    end

    c1 = mean(bnd);
    G1 = diag((bnd(:,1)-bnd(:,2))./2);
    c = [c1; c1];
    G = blkdiag(G1,G1);
    A = [G1, -G1];
    b = shift;%<== difference between generator*factor = shift
    vset = true(1,size(G,2));
    keys.dims = dimKeys;
    keys.factors = memZono.genKeys(lbl,1:size(G,2));
    keys.cons = memZono.genKeys(lbl,1:size(A,1));
    F_shift = memZono(G,c,A,b,vset,keys);
end


function y = FmEval(Fm,x,inDims,outDims)%<= dims = {x1,x2,y}
    arguments
        Fm
        x
        inDims = 'x';
        outDims = 'y';
    end
    X = memZono(x,inDims);
    Y_Fm = Fm.and(X,'eval');
    y = Y_Fm
end