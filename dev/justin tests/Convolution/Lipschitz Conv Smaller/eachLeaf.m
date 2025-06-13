% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v] = eachLeaf(obj,optSolver)

% Determine number of non-empty constrained zonotopes
[leaves] = getLeaves(obj,optSolver);
if isempty(leaves)
    warning('zonoLAB:EmptyZonotope','Hybrid zonotope is empty and cannot be plotted.')
    v = [];
    return
end
nLeaves = size(leaves,2);

optPlot = plotOptions('Display','off','SolverOpts',optSolver);
v = [];
nVerts = zeros(nLeaves,1);
waitbarHandle = waitbar(0,['Plotting hybrid zonotope with ',num2str(nLeaves),' leaves.']);
%  If \xib denotes the i^th column of leaves, then the corresponding
%  non-empty constrained zonotope is of the form:
%  Z = { (c + Gb \xib) + Gc \xic | ||\xic||_inf <= 1, Ac \xi = b - Ab \xib }
emptyLeaves = false;
for i = 1:nLeaves
    Zi = conZono(obj.Gc,obj.c+obj.Gb*leaves(:,i),obj.Ac,obj.b-obj.Ab*leaves(:,i));
    [vi] = findVF(Zi,optPlot.SolverOpts);
    v = [v; vi];
    waitbar(i/nLeaves,waitbarHandle)
end
close(waitbarHandle)

% Check for unplotted conZono leaves
if emptyLeaves
    warning('zonoLAB:Tolerance','Some leaves of the hybrid zonotope did not plot and may be caused by constraints that are nearly redundant (close to the plotting/optimization tolerances). Check the validity of other leaves.')
end

end