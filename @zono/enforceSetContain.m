% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Enforces constraints to check if the constrained zonotope X is contained within the constrained zonotope Y.
%   Syntax:
%       [prob] = enforceSetContain(X, Y, prob, varargin)
%   Inputs:
%       X - zonotope
%       Y - zonotope
%       prob - MATLAB optimproblem
%       varargin - string or character input to add as a suffix for variables if doing nested optimization to
%       enforce uniqueness (default is '1'.)
%
%       
%   Outputs:
%       prob - returns the optimproblem containing the constraints necessary to check for constrained zonotope set containment.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [prob] = enforceSetContain(X,Y,prob,suffix)
    arguments
        X zono
        Y zono
        prob
        suffix string = '1'
    end

% %check for user input of suffix variable, otherwise 1
% if length(varargin) > 0
%     suffix = string(varargin);
% else
%     suffix = '1';
% end

% sizes
n_x = size(X.G, 2);
n_y = size(Y.G, 2);

% Create optimization variables
gamma = optimvar(sprintf('gamma_%s',suffix), n_y, n_x);
beta = optimvar(sprintf('beta_%s',suffix), n_y);
g = optimvar(sprintf('g_%s',suffix), n_y, n_x);
b = optimvar(sprintf('b_%s',suffix), n_y);

% constraints
if isempty(prob.Constraints), prob.Constraints = struct; end

prob.Constraints.(sprintf('gamma_%s',suffix)) = X.G == Y.G*gamma;
prob.Constraints.(sprintf('beta_%s',suffix)) = Y.c - X.c == Y.G*beta;
prob.Constraints.(sprintf('gamma_abs_1_%s',suffix)) = g >= gamma;
prob.Constraints.(sprintf('gamma_abs_2_%s',suffix)) = g >= -gamma;
prob.Constraints.(sprintf('beta_abs_1_%s',suffix)) = b >= beta;
prob.Constraints.(sprintf('beta_abs_2_%s',suffix)) = b >= -beta;
prob.Constraints.(sprintf('gamma_beta_%s',suffix)) = g*ones(n_x,1) + b <= ones(n_y,1);

end

