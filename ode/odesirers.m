% This plots the error space of the function within the pre-specified 
% values of beta and gamma. Uses odeBDSIR function to calculate the error.

% Resolution of the error space (1000 by 1000)
RES = 4;

% Limit of beta and gamma
PAR_LIMIT = 1;

% Top number of parameters to obtain
TOP = 10;

% Generate beta and gamma in the error space
beta = linspace(0, PAR_LIMIT, RES);
gamma = beta;

% Empty vector for error and rnaught
error = zeros(RES);
rnaught = zeros(RES);

% Calculate error and rnaught
for i = 1:RES
    for j = 1:RES
        par.beta = beta(i);
        par.gamma = gamma(j);
        error(i,j)= ode_sir(par).fitMeasure;
        rnaught(i,j) = par.beta / par.gamma;
    end
end

% Plots the error surface
error = log(error);
surf(beta, gamma, error)
colorbar
caxis([0 1000])
shading flat
title('Error Space of Deterministic SIR model')
xlabel('beta')
ylabel('gamma')
zlabel('log(error)')

% Obtain the top 10 parameters with minimum error
[errorSort, errorInd] = sort(error(:), 1, 'ascend');        %sorts errors from lowest to highest
min_errors = errorSort(1:TOP);                              %obtains the 10 lowest errors
[ind_row, ind_col] = ind2sub(size(error), errorInd(1:TOP)); %fetch indices of the 10 lowest errors

% Stores best parameters and the corresponding r-naught and error to one matrix
min_param = zeros(TOP, 4);
for i = 1:TOP
        min_param(i,1) = beta(ind_row(i));
        min_param(i,2) = gamma(ind_col(i));
        min_param(i,3) = min_param(i,1) / min_param(i,2);
        min_param(i,4) = min_errors(i);
end

% Plot top 10 beta
boxplot(min_param(:,1))
title('Beta of the top 10 Parameter Estimates with Minimum Error')
ylabel('Beta')

% Plot top 10 gamma
boxplot(min_param(:,2))
title('Gamma of the top 10 Parameter Estimates with Minimum Error')
ylabel('Gamma')

% Plot top 10 r-naught
boxplot(min_param(:,3))
title('R-naught of the top 10 Parameter Estimates with Minimum Error')
ylabel('R0')


function out = ode_sir(par)
% ode_sir calculates the SIR projections using the given parameters

    % Read the NCR reported case data
    dat = readtable("ncr-reported.xlsx");

    % Initial conditions (in terms of proportion)
    popN_raw = 7e6;
    popN = 1;                               
    I0 = mean(dat.reported(1:5))/popN_raw;       
    x0 = [popN-I0 I0 0];

    % Define the timespan (in days)
    tspan = caldays(between(dat.date(1), dat.date(end), 'Days')); 
    timeDomain = [0 tspan];

    % Simulate the SIR projections
    options = odeset( ...
        'RelTol',1e-10, ...     % High accuracy
        'AbsTol',1e-10 ...      % Consider solutions very close to zero
        );
    [out.T, out.Y] = ode15s(@model, timeDomain, x0, options);

    % Calculates SIR in terms of raw values (not proportion)
    out.Y = popN_raw*out.Y;
    out.Trep = caldays(between(dat.date(1), dat.date, 'Days'));
    out.Yest = interp1(out.T, out.Y(:,2), out.Trep);
    out.Sus = interp1(out.T, out.Y(:,1), out.Trep);
    out.Rec = interp1(out.T, out.Y(:,3), out.Trep);
    out.Yobs = dat.reported;

    % Calculates the error value
    out.fitMeasure = sqrt(mean((out.Yobs - out.Yest).^2));
    
    function dx = model(~,x)
    % model defines the SIR model

        % Input state values
        popS = x(1);
        popI = x(2); 

        % Declare equation terms
        infection = par.beta*popI*popS; 
        removal = par.gamma*popI;
        dx = zeros(3,1);
        dx(1) = -infection;
        dx(2) = infection - removal;
        dx(3) = removal;
    end
end