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