function out = odesir_model(par, graph)
% ode_sir calculates the SIR projections using the given parameters

    % Number of days to project
    if graph == 0 % No projection yet for data fitting process
        DAYS_PROJECT = 0;
    elseif graph == 1
        DAYS_PROJECT = 300;
    end

    % Read the NCR reported case data
    bdsir = readtable("ncr-bdsir-infected.csv");

    % Interpolate values
    dat.date = datetime(2020,03,02,0,0,0):datetime(2020,7,26,0,0,0);
    dat.bdsir = interp1(bdsir.bdsir_date', bdsir.bdsir_infected', dat.date);

    % Initial conditions (in terms of proportion)
    popN_raw = 7e6;
    popN = 1;                               
    I0 = mean(dat.bdsir(1:5))/popN_raw;       
    x0 = [popN-I0 I0 0];

    % Define the timespan (in days)
    tspan = caldays(between(dat.date(1), dat.date(end), 'Days')) + DAYS_PROJECT; 
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
    if graph == 1 
        out.Trep_obs = out.Trep;  % Create time axis for observations
        out.Trep = [out.Trep  tspan-DAYS_PROJECT+1:tspan-DAYS_PROJECT+300]; % Extend simulation for projection
    end
    out.Yest = interp1(out.T, out.Y(:,2), out.Trep);
    out.Sus = interp1(out.T, out.Y(:,1), out.Trep);
    out.Rec = interp1(out.T, out.Y(:,3), out.Trep);
    out.Yobs = dat.bdsir;

    % Calculates the error value
    if graph == 1 % Do not include projected values in calculating error for graphing
        out.fitMeasure = sqrt(mean((out.Yobs - interp1(out.T, out.Y(:,2), caldays(between(dat.date(1), dat.date, 'Days')))).^2));
    elseif graph == 0
        out.fitMeasure = sqrt(mean((out.Yobs - out.Yest).^2));
    end
    
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