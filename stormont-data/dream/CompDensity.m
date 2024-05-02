function [p,log_p] = CompDensity(x,MCMCPar,Measurement,ModelName,Extra,option)
% This function computes the density of each x value

p = []; log_p = [];

% Loop over the individual parameter combinations of x
for ii = 1:size(x,1),
    % Call model to generate simulated data
    evalstr = ['ModPred = ',ModelName,'(x(ii,:),Extra);']; eval(evalstr);
    
    if option == 1, % Model directly computes posterior density
        p(ii,1) = ModPred; log_p(ii,1) = log(p(ii,1));
    end;

    if option == 2, % Standard log-likelihood function
        % Calculate the error residual
        Err = ( Measurement.MeasData(:) - ModPred(:) ); 
        % Derive the log likelihood
        if size(Measurement.Sigma,1) == 1,
            log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - ( Measurement.N / 2) * log( Measurement.Sigma ) - 1/2 * Measurement.Sigma^(-2) * sum ( Err.^2 );
        else
            log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - ( sum ( log( Measurement.Sigma ) ) / 2) - 1/2 * sum ( ( Err./Measurement.Sigma ).^2);
        end;
        % And retain in memory
        p(ii,1) = log_p(ii,1);
    end;

    if option == 3, % Model computes output simulation
        %Err = (Measurement.MeasData(:)-ModPred(:));
        Err = (Measurement.MeasData(:)-ModPred(:)).*Measurement.Weights;
        % Derive the sum of squared error
        SSR = sum(Err.^2);
        % And retain in memory
        p(ii,1) = -SSR;
	log_p(ii,1) = -0.5 * SSR;
    end;

    if option == 4, % Model directly computes log posterior density
        p(ii,1) = ModPred; log_p(ii,1) = p(ii,1);
    end;

    if option == 5, % Similar as 3, but now weights with the Measurement Sigma
        % Defime the error
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the sum of squared error
        SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
        % And retain in memory
        p(ii,1) = -SSR; log_p(ii,1) = -0.5 * SSR;
    end;

    if option == 8, % Generalized log likelihood (GL)
        % Extract statistical model parameters
        par = Extra.fpar;               % fixed parameters
        par(Extra.idx_vpar) = x(ii,:);  % variable parameters
        par = par';                     % make it a column vector
        statpar = par(end-10:end);
        % Compute the log-likelihood
        log_p(ii,1) = GL('est',statpar,ModPred,Measurement.MeasData);
        % And retain in memory
        p(ii,1:1) = log_p(ii,1);
    end;

    if option == 9, % Whittle likelihood function
        
        % Calculate the log-likelihood using spectral density 
        [log_L] = Whittle_logL(Measurement,ModPred);        
        
        % Now store in memory
        p(ii,1) = log_L; log_p(ii,1) = p(ii,1);
                
    end;

end;
