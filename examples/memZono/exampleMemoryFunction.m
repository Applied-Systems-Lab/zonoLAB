%% Example Function Creation
clear; 
% close all;

%% Function Setup
n = 11;
x0 = 10;
bd = [-x0,x0];
f = @(x) sin(x);
% f = @(x) atan(x);

f_SOS = makeSOSfunc(f,n,bd);
F_SOS = memZono(f_SOS,{'x','y'},'sos');

%% Account for error
errSOS = makeSOSerror(f,n,bd);
errSOS_max = max(abs(errSOS),[],'all');
f_Err = zono([0;errSOS_max],zeros(2,1));
F_Err = memZono(f_Err,{'x','y'},'err');

%% Construct function
F = f_SOS+f_Err;
Fm = F_SOS + F_Err;


if false
%% plotting
figure;
t = tiledlayout("vertical");
% f(x)
dx = 0.1;
x = bd(1):dx:bd(2);
y = f(x);
% sos
x_sos = linspace(bd(1),bd(2),2*n+1);
y_sos = f(x_sos);

%% Plot full region
ax = nexttile; hold on;
Fm.plot({'x','y'},'b',0.5);
plot(x_sos,y_sos, 'k--','DisplayName','f_{SOS}(x)');
plot(x,y,'r--','DisplayName','f(x)');
legend


%% X intersections
% dx_sos = 0.25;%diff(x_sos([1 2]));
% for i = 1:numel(x_sos)
%     X_{i} = zono([dx_sos/4],x_sos(i));
% end
% Xall = [X_{:}];
% X = X_{1};
% for i = [3,5,7,9,11]
%     X = union(X,X_{i});
% end
% % for i=; X = union(X,Xall(i)); end

X = zono([pi/2],0);
Xm = memZono(X,'x','X');

FXm = and(Fm,Xm,'eval');
ax = nexttile; hold on;
plot(FXm,{'x','y'},'r',0.5);
plot(x,y,'r--');

% FXm2 = and(F,X,[1,0])
% plot(FXm2,'r',0.1)



% Y intersections
Y = zono(0.25,0.25);
Ym = memZono(Y,'y','Y');

FYm = and(Fm,Ym,'eval');
ax = nexttile; hold on;
plot(FYm,{'x','y'},'b',1);
plot(x,y,'r--');

%
Ym2 = memZono(union(Y+0.25,-1*(Y+0.25)),'y','Y2');
FYm2 = and(Fm,Ym2,'eval');
ax = nexttile; hold on;
plot(FYm2,{'x','y'},'b',1);
plot(x,y,'r--');

end



%% Shifting
figure; 
t = tiledlayout('vertical');
ax = nexttile; hold on;

%% Domain Shift
xshift = pi/2;
shift_X = makeFshift(xshift,2*bd);
shift_Xm = memZono(shift_X,{'in','out'},'xshift');
Fm_shiftX = and(...
    relabelDims(shift_Xm,{'in','out'},{'x_xshift','inter'}),...
    relabelDims(Fm,{'x','y'},{'inter','f_xshift'}),...
'shift');

%% Range shift
yshift = 0.5;
shift_Y = makeFshift(yshift,[-2 2]);
shift_Ym = memZono(shift_Y,{'in','out'},'yshift');
Fm_shiftY = and(...
    relabelDims(shift_Ym,{'in','out'},{'inter','f_yshift'}),...
    relabelDims(Fm,{'x','y'},{'x_yshift','inter'}),...
'shift');

%% Ploting
% Fall = [Fm; Fm_shiftX; Fm_shiftY];
plot(Fm,{'x','y'},'k','0.5');
plot(Fm_shiftX,{'x_xshift','f_xshift'},'r',0.5);
plot(Fm_shiftY,{'x_yshift','f_yshift'},'b',0.5);


%% adding F and shift F
ax = nexttile; hold on;
plot(Fm,{'x','y'},'k','0.5');
plot(Fm_shiftX,{'x_xshift','f_xshift'},'r',0.5);
Fm_plusShift = Fm + relabelDims(Fm_shiftX,{'x_xshift','f_xshift'},{'x','y'});
Fm_plusShift.plot({'x','y'},'g',0.5);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeSOSfunc
function f_SOS = makeSOSfunc(f,n,bd)
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
    
    f_SOS = hybZono(Gc,Gb,c,Ac,Ab,b);
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
function F_shift = makeFshift(shift,bnd)
    arguments
        shift
        bnd = [-10,10]
    end

    c1 = mean(bnd);
    G1 = diag((bnd(:,1)-bnd(:,2))./2);
    c = [c1; c1];
    G = blkdiag(G1,G1);
    A = [G1, -G1];
    b = shift;%<== difference between generator*factor = shift
    F_shift = conZono(G,c,A,b);
end