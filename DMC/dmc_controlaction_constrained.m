function du = dmc_controlaction_constrained(dmc, yr, dup, yp, up, ulimits, dulimits, ylimits)
%Given a Dynamic Matrix Control (DMC) structure and control inputs
%calculates the control action which minizes DMC's cost function.
%
% Inputs:
%         dmc - a DMC structure setup by dmc_setup function
%         yr - control reference
%         dimension: Ny-N1+1 x output_count
%         dup - set of past delta (variation of) process inputs
%         dimension: N-N1 x input_count
%         yp - set of past process outputs
%         dimension: output_count
%         up - set of past (t-1) process inputs, necessary for limits
%         dimension: input_count
%         ulimits - set of lower and upper limits on u on the form
%         [u1min, u1max; u2min, u2max; ... ]
%         dimension: input_count x 2
%         dulimits - set of lower and upper limits on delta_u on the form
%         [du1min, du1max; du2min, du2max; ... ]
%         dimension: input_count x 2
%         ylimits - set of lower and upper limits on y on the form
%         [y1min, y1max; y2min, y2max; ... ]
%         dimension: output_count x 2
%

% Calculate free response and error estimation
yp = repmat(yp, dmc.Ny-dmc.N1+1, 1);

f = dmc.dG*dup(:) + yp(:);
e = yr(:) - f;

% Check inputs
if(~exist('ulimits', 'var'))
    ulimits = [];
elseif(~isempty(ulimits))
    if(isempty(up))
        error('DMC: Input constraints requires past u(k-1) inputs');
    end
    
    %  du(t) <= max - u(t-1)
    %  du(t) >= min - u(t-1)
    tmp = repmat({[tril(ones(dmc.Nu));-tril(ones(dmc.Nu))]},dmc.input_count,1);
    Au = blkdiag(tmp{:}); % square matrix (Nu * input_count)
    
    up = up(:);
    ulimits = ulimits - repmat(up, 1, dmc.input_count);
    bu = [repmat(ulimits(:,2), 1, dmc.Nu)'; -repmat(ulimits(:,1), 1, dmc.Nu)'];
    bu = bu(:);
end

if(~exist('dulimits', 'var'))
    dulimits = [];
elseif(~isempty(dulimits))
    %  du(t) <= max
    %  du(t) >= min
    tmp = repmat({[tril(ones(dmc.Nu));-tril(ones(dmc.Nu))]},dmc.input_count,1);
    Adu = blkdiag(tmp{:}); % square matrix (Nu * input_count)
    
    bdu = [repmat(dulimits(:,2), 1, dmc.Nu)'; -repmat(dulimits(:,1), 1, dmc.Nu)'];
    bdu = bdu(:);
end

if(~exist('ylimits', 'var'))
    ylimits = [];
elseif(~isempty(ylimits))
    %  du(t) <= max
    %  du(t) >= min
    Ay = [dmc.G; -dmc.G]; % square matrix (Nu * input_count)
   
    by = [repmat(ylimits(:,2), 1, dmc.Ny-dmc.N1+1)'; -repmat(ylimits(:,1), 1, dmc.Ny-dmc.N1+1)'];
    by = by(:) - [f; -f];
end

% Calculate control action increments
if(isempty(ulimits) && isempty(dulimits) && isempty(ylimits))
    du = dmc.K*e;
else
    % Assemble A/b restrictions
    A = [];
    b = [];
    if(~isempty(ulimits))
        A = [A; Au]; b = [b; bu];
    end
    
    if(~isempty(dulimits))
        A = [A; Adu]; b = [b; bdu];
    end
    
    if(~isempty(ylimits))
        A = [A; Ay]; b = [b; by];
    end
    
    % Use quadratic programming to solve the optimization problem
    opts = optimset('Display','off');
    [du,~,~,~,~] = quadprog(dmc.H,-dmc.G'*e,A,b,[],[],[],[],[],opts);
    
    % Take only the first du(t) of each input
    du = du(1:dmc.Nu:end,1);
end

end

