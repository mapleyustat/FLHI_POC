function du = dmc_controlaction(dmc, yr, dup, yp)
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
%

% Calculate free response and error estimation
yp = repmat(yp, dmc.Ny-dmc.N1+1, 1);

f = dmc.dG*dup(:) + yp(:);
e = yr(:) - f;

% Calculate control action increments
du = dmc.K*e;

end

