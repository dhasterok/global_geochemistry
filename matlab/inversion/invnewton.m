function [m,rms,ci95] = invnewton(func,x,d,m0,varargin)
% INVNEWTON - Newtons inversion method
%
%   [M,RMS,CI95] = INVNEWTON(FOO,X,D,M0) computes the model, M, with
%   constants used to fit the data, D, with independent variable(s) in X.
%   Each independent variable associated with D should be a single column
%   of LENGTH(D).  The forward model and Frechet matrix computations are
%   handeled by the function named in FOO.  See below for an example.
%   The initial guess, M0, is a column vector.  The RMS misfit and 95%
%   confidence intervals (CI95) may also be returned.
%
%   [M,RMS] = INVNEWTON(FOO,X,D,M0,wd,wx) adds the effect of data (wd) and
%   model (wx) weighting to find the "best- fitting" model.  Adding model
%   weights requires an additional function to be added within FOO.
%
%   If necessary, a tolerance, TOL, and maximum number of iterations,
%   MAXIT, can be added.
%
%   [M,RMS] = INVNEWTON(FOO,X,D,M0,wd,wx,sig_m,TOL,MAXIT)
%
%   In order for INVNEWTON to work, a function named in FOO must be
%   included in the directory.  FOO should be called by [DM,F] = FOO(X,M0)
%   where DM is the forward model and F is the Frechet matrix.  The
%   function should be stored in the current directory.
%
%   Example function FOO = 'model':
% 
%   File: model.m
%   
%   function [dm,F] = model(m,x)
%   
%   dm = forward(m,x);
%   F = frechet(m,x);
%   
%   return
%   
%   function dm = forward(m,x)
%       dm = m(1)*(298./x).^m(2);
%   return
%   
%   function F = frechet(m,x)
%       F = [(298./x).^m(2) ...
%            m(1)*(298./x).^m(2).*log(298./x)];
%   return
%
%   Note: computing the Frechet matrix can be very costly.  If you wish to
%   call only the forward model to compute a model solution, it might be
%   wise to use VARARGOUT and NARGOUT to optionally compute the Frechet
%   matrix only when two variables are returned.
%
%   Example:
%
%       x = [773:100:1273]';
%       d = [2.93 2.72 2.66 2.54 2.46 2.34]';
%       m0 = [4; 0.5];
%
%       [m,rms] = invnewton('model',x,d,m0);
%

% Original: 13-Feb 2009 by D. Hasterok
%
% Revised: 28-Jul 2010 by D. Hasterok.  Changed exit on increase in RMS to
% a reduction in the step size, reducing the chances of steping beyond the
% minima that would otherwise result in an increase in the RMS.  The step
% size reduction is limited to 8 before accepting defeat.
%
% Revised: 7-Oct 2014 by D. Hasterok.  Added model weights (requires an
% additional function to be added to the model.m file.

m0 = m0(:);
d = d(:);
if min(size(x) == 1)
    x = x(:);
end
% Nd number of data
% Np number of parameters
% Nm number of model parameters
Nm = length(m0);
[Nd Np] = size(x);

if length(d) ~= Nd
    error('(INVNEWTON): number of data do not match number of parameter pairs.');
end

% End conditions
maxit = 1e2;
tol = 0.5*(max(d) - min(d))*1e-6;
%tol = 1e-4;

% Data weights
W = 1;

weights = 0;
if nargin > 3 | nargin < 9
    if nargin > 4
        if nargin == 5
            error('(INVNEWTON) Incorrect number of input arguments.');
        end

        % Weighting
        sigma_d = varargin{1};
        sigma_x = varargin{2};
        if ~isempty(sigma_x) | ~isempty(sigma_d)
            W = eye(Nd);

            % Were data weights supplied?
            if length(sigma_d) ~= 1 & length(sigma_d) ~= Nd & ~isempty(sigma_d)
                error('(INVNEWTON): length(sigma_d) == number of data.');
            elseif length(sigma_d) == 1
                sigma_d = sigma_d*ones([Nd 1]);
            end

            % Were parameter weights supplied?
            [nr,nc] = size(sigma_x);
            if nr ~= Nd & nr ~= 1 & nc ~= Np & ~isempty(sigma_x)
                error('(INVNEWTON): size(sigma_x) == size(x).');
            elseif nr == 1
                temp = sigma_x;
                sigma_x = zeros([Nd Np]);
                for i = 1:Np
                    sigma_x(:,i) = temp(i);
                end
            end

            if ~isempty(sigma_d)
                weights = 1;

                % compute variance of input data from standard deviations
                var_d = sigma_d(:).^2;

                % if there are no parameter weights, then the weighting matrix
                % is diagonal with the elements as reciprocal of variance
                if isempty(sigma_x)
                    W = diag(1./var_d, 0);
                end
            else
                var_d = 0;
            end

            if ~isempty(sigma_x)
                weights = 2;

                % compute variance of input parameters from standard deviations
                var_x = sigma_x.^2;
            end
        end

        % additional arguements for tolerance and maximum number of iterations
        if nargin > 6
            if nargin == 7
                error('(INVNEWTON) Incorrect number of input arguments.');
            end
            % Was tolerance and/or maximum number of iterations supplied?
            tol = varargin{2};
            if nargin == 7
                maxit = varargin{3};
            end
        end
    end
else
    error('(INVNEWTON) Incorrect number of input arguments.');
end

% Number of samples to fit
if weights == 0
    N = Nd;
elseif weights == 1
    N = sum(diag(W));
end

iter = 1;
n = 1;
while 1
    if iter > maxit
        warning('  (INVNEWTON) Number of iterations exceeds allowable.');
        break
    end

    % Compute starting solution and set initial RMS
    if iter == 1
        if weights == 2
            % if there are uncertainties associated with the input parameters,
            % then the weighting matrix is determined by the reciprocal of
            % propagation of errors (this must be updated every step)

            % compute derivatives of forward w.r.t. input parameters
            [dm,F,Kx] = feval(func,m0,x);

            Md = diag(var_d,0);
            Mx = diag( sum(Kx .* var_x, 2) , 0);

            % weighting matrix
            W = inv( Md + Mx );
            N = sum(diag(W));
        else
            % if the weighting matrix is already determined
            [dm,F] = feval(func,m0,x);
        end
        rmsold = sqrt(sum(W*(dm - d).^2)/N);
        fprintf('  Starting Model: %3i   RMS: %.4f\n',iter-1,rmsold);
    end

    % Compute updated model (unregularized)
    warning off;
    m = m0 - n*inv(F'*W*F)*F'*W*(dm - d);
    warning on;

    % Compute forward solution with new model
    % and Frechet matrix for next iteration
    [dm,F] = feval(func,m,x);
    if ~isempty(find(isnan(dm))) | ~isempty(find(isnan(F)))
        m(:) = m0;
        rms = rmsold;
        %ci95 = 1.96*sqrt(diag(inv(F'*F)));
        warning('  WARNING (model): Found nans in dm.');
        break
    end

    rms = sqrt(sum(W*(dm - d).^2)/N);
    fprintf('  Iteration:      %3i   RMS: %.4f\n',iter,rms);

    % Check termination conditions
    if abs(rms - rmsold) < tol
        break
    end

    % Update values and counter
    if rms < rmsold
        % step-size was appropriate, move on to the next step.
        n = 1;
        m0 = m;
        rmsold = rms;
        iter = iter + 1;
    else % rms >= rmsold
        % assume the step-size was too large, reduce the model step and
        % try again.
        n = 0.5*n;
        if n < 0.0078125
            warning('  WARNING (invnewton): RMS increased quitting...');
            break
        end
    end
end

ci95 = 1.96*sqrt(diag(inv(F'*F)));

return