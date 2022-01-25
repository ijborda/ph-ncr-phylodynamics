% This script find the beta and gamma that best fits the NCR reported 
% cases using the classical SIR model

% Define random seed
rng(3453);

% Number of simulation
simNumber = 1;

% Parameter limits
BETA_LIM = 4;
GAMMA_LIM = 4;

for i = 1:simNumber
  % Initialize random initial beta and gamma
  initialBeta = BETA_LIM * rand; 
  initialGamma = GAMMA_LIM * rand;
 
  % Find best fitting gamma and beta
  options = optimset ('TolFun', 1e-50, 'TolX', 1e-50, 'MaxFunEvals', 1e100);
  [bestPar, bestError] = fminsearch(@ode_sir_error, [initialBeta initialGamma], options);
  
  % Store simulation result
  sim.initialBeta(i) = initialBeta;   
  sim.initialGamma(i) = initialGamma;   
  sim.beta(i) = bestPar(1);       
  sim.gamma(i) = bestPar(2);           
  sim.r(i) = bestPar(1) / bestPar(2); 
  sim.error(i) = bestError;  
end


% Plot best fitting beta and gamma
scatter(sim.beta, sim.gamma, 100, '.');
xlim([0 BETA_LIM]);
ylim([0 GAMMA_LIM]);
title(sprintf('Best fitting Beta and Gamma (simNumber = %i)', simNumber));
xlabel('Beta');
ylabel('Gamma');
xtickformat('%,.2f');
ytickformat('%,.2f');

figure();

% Boxplot for best r
boxplot(sim.r);
title(sprintf('Best fitting Basic Reproductive Numbers (simNumber = %i)', simNumber));
ytickformat('%,.2f');
ylim([0 BETA_LIM]);

figure();

% Graph infected trajectory against NCR cases using parameters with least error 
[errorSort, errorInd] = sort(reshape(sim.error,[],1), 1, 'ascend'); 
par.beta = sim.beta(errorInd(1));
par.gamma = sim.gamma(errorInd(1));
out = odesir_model(par, 1);
plot(out.T, out.Y(:,2), out.Trep, out.Yest, 'o', out.Trep_obs, out.Yobs);
legend({'Model (model)','Infected (projected)','Infected (reported)'},'Location','best');
title(sprintf('Projected vs. Reported \n (Beta = %.4f, Gamma = %.4f, R0 = %.2f, Error = %.2f)', par.beta, par.gamma, par.beta / par.gamma, sim.error(errorInd(1))));

function error = ode_sir_error(p)
% Returns error value in the SIR projection given beta and gamma

    % Returns very large value if beta or gamma is negative     
    if sum(p < 0) > 0
        error = 1e6;

    % Otherwise, calculates error using ode_sir function 
    else
        par.beta = p(1);
        par.gamma = p(2);
        error = odesir_model(par, 0).fitMeasure;
    end
end