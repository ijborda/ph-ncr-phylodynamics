% Graph infected trajectory against NCR cases using parameters with least error 
par.beta = 0.0891;
par.gamma = 0.0611;
out = odesir_model_projection(par, 1);
plot(out.T, out.Y(:,2), out.Trep, out.Yest, 'o', out.Trep_obs, out.Yobs);
legend({'Model (model)','Infected (projected)','Infected (reported)'},'Location','best');
title(sprintf('Projected vs. Reported \n (Beta = %.4f, Gamma = %.4f, R0 = %.2f)', par.beta, par.gamma, par.beta / par.gamma));