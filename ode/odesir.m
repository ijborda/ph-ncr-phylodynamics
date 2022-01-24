% This plots the error space of the function within the pre-specified 
% values of beta and gamma. Uses odeBDSIR function to calculate the error.

% Resolution of the error space (1000 by 1000)
RES = 5;

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
iter = 1;
for i = 1:RES
    for j = 1:RES
        fprintf('Iteration: %d\n', iter);    %for monitoring of progress
        par.beta = beta(i);
        par.gamma = gamma(j);
        error(i,j)= odesir_model(par).fitMeasure;
        rnaught(i,j) = par.beta / par.gamma;
        iter = iter + 1;
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

figure ();

% Plot top 10 gamma
boxplot(min_param(:,2))
title('Gamma of the top 10 Parameter Estimates with Minimum Error')
ylabel('Gamma')

figure ();

% Plot top 10 r-naught
boxplot(min_param(:,3))
title('R-naught of the top 10 Parameter Estimates with Minimum Error')
ylabel('R0')

figure ();

% Graph infected trajectory against NCR cases using parameters with least error 
par.beta = min_param(1,1);
par.gamma = min_param(1,2);
out = odesir_model(par);
plot(out.T, out.Y(:,2), out.Trep, out.Yest, 'o', out.Trep, out.Yobs);
legend({'Model (model)','Infected (projected)','Infected (reported)'},'Location','best');
title(sprintf('Projected vs. Reported \n (Beta = %.4f, Gamma = %.4f, R0 = %.2f, Error = %.2f)', par.beta, par.gamma, par.beta / par.gamma, min_param(1,4)));