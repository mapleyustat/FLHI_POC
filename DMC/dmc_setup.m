function [dmc] = dmc_setup(N, N1, Ny, Nu, lambda, step_model, ulimits, dulimits, ylimits)
%% Setups Dynamic Matrix Control (DMC) main variables
%
% Inputs:
%         N - step model horizon
%         N1 - start of prediction horizon
%         Ny - prediction horizon
%         Nu - control horizon
%         lambda - control weighing factor
%

samples = size(step_model, 1);
input_count = size(step_model, 2);
output_count = size(step_model, 3);

if(samples ~= N)
    error('Step model sample count does not match step model horizon');
end

dmc.G = zeros((Ny-N1+1)*output_count, Nu*input_count);
dmc.dG = zeros((Ny-N1+1)*output_count,(N-N1)*input_count);
dmc.K = zeros(input_count,(Ny-N1+1)*output_count);

%% Create and fill the dynamic matrix (G)
for outputIndex = 1:output_count
    for inputIndex = 1:input_count
        % Get the step coefficients for this transfer function
        g = step_model(:, inputIndex, outputIndex);
        
        % (Ny-N1+1)xNu sized
        g_idx = repmat([N1:Ny]', 1, Nu) - repmat([0:Nu-1], Ny-N1+1, 1);
        g_idx_null = g_idx < 1;
        g_idx(g_idx < 1) = 1;
        G = g(g_idx);
        G(g_idx_null) = 0;

        % Create and fill the deltag matrix
        % (Ny-N1+1)x(N-N1) sized
        deltag_fixed_idx = repmat(N1:N-1, Ny-N1+1, 1);
        deltag_for_idx = repmat(N1+1:N+1-1, Ny-N1+1, 1) + repmat([0:Ny-N1]', 1, N-N1);
        % if g(j+1) is beyond the horizon of the step response, it's considered asymptotically stable
        deltag_for_idx(deltag_for_idx > N) = N;

        dG = g(deltag_for_idx) - g(deltag_fixed_idx);
        
        dmc.G(1+(outputIndex-1)*(Ny-N1+1):outputIndex*(Ny-N1+1), 1+(inputIndex-1)*Nu:inputIndex*Nu) = G;
        dmc.dG(1+(outputIndex-1)*(Ny-N1+1):outputIndex*(Ny-N1+1), 1+(inputIndex-1)*(N-N1):inputIndex*(N-N1)) = dG;
    end
end


%% Calculate the gain vector K from [G'*G + lambda*I]^(-1)*G'
% Calculate the inverse portion of the gain
dmc.H = dmc.G'*dmc.G + lambda*eye(Nu*input_count);

% Calculate a gain matrix K for each output
% first it's Nu*input_count by (Ny-N1+1)*output_count
% then input_count by (Ny-N1+1)*output_count
dmc.Kfull = dmc.H\dmc.G';
dmc.K = dmc.Kfull(1:Nu:end, :); % use only the first du(t) of each input

%% Finish setting up the structure
dmc.N = N;
dmc.N1 = N1;
dmc.Ny = Ny;
dmc.Nu = Nu;
dmc.lambda = lambda;

dmc.input_count = input_count;
dmc.output_count = output_count;

end
