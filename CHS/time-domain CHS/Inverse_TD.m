function [xfit,xfit_s,datafit,chi2r] = Inverse_TD(dynx,ctHb,data,data_err,Fs)
    % Inverse model for time-domain coherent hemodynamics spectroscopy
    % Inputs:
    %       dynx - CHS dynamic parameters: blood volume contributions from
    %               the arteries, capillaries, and veins to the total blood volume
    %              dynx = [CBVa/CBV, CBVc/CBV, CBVv/CBV].
    %       ctHb -  hemoglobin concentration in blood. Unit is microM/unit of blood volume
    %       data   -   Temporal trace of measured absolute oxyhemoglobin and 
    %                  deoxyhemoglobin concentrations [Ostep,Dstep]. 
    %               Units are microM.
    %       data_err  - Error of measured O and D. Units are microM.
    %       Fs  -   Sampling rate. Unit is Hz
    % Outputs:
    %       xfit  -   CHS baseline parameters
    %                 x = [tc tv F(c)CBV0c/CBV0 fc k]. Units are [s s - Hz -].  
    %       xfit_s -  estimated errors of fitted CHS baseline parameters
    %                   obtained by using bootstrapping method
    %       datafit - Temporal trace of absolute oxyhemoglobin concentration
    %               obtained from the fit. Units are microM.
    %       chi2r   - reduced chi-squared of the fit.

    % Written By Jana Kainerstorfer
    % Modified by Kristen Tgavalekos
    % Last modified by Thao Pham (February 2020)

    %% Setups 
    Ostep = data(:,1);
    Dstep = data(:,2);
    Tstep = Ostep+Dstep;
    T0 = Tstep(1);
    
    %% Fit
    % [tc tv FCBV0c/CBV0 fc k]
    x0 = [1  5   0.3  0.03   2];  % Initial parameters 
    lb = [0  0   0    0      1];  % lower bound of the fit
    ub = [2  10  0.7  0.2    10]; % upper bound of the fit
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    
    fwd_func = @(x)forward_TD(x,dynx,Tstep,T0,ctHb,Fs);
    [xfit,chi2r] = fmincon(@(x0)costfunc(fwd_func,x0,data,data_err),x0,[],[],[],[],lb,ub,[],options); 
    [Ofit,Dfit] = feval(fwd_func,xfit);
    datafit = [Ofit,Dfit]; 
    
    %% Calculate errors in the fitted parameters
    nboot = 50;
    rng(14,'twister');
    
    % Compute residual
    residuals = data-datafit;
    [~,bootIndices] = bootstrp(nboot,[],residuals);
    
    % Start bootstrapping
    xfitboot = zeros(nboot,length(xfit));
    fprintf('Bootstrapping:');
    for i = 1:nboot
        fprintf(' %i,',i);
        bootvals_O = residuals(bootIndices(:,i),1);
        bootvals_D = residuals(bootIndices(:,i),2);
        databoot(:,1) = datafit(:,1)+bootvals_O;
        databoot(:,2) = datafit(:,2)+bootvals_D;
        Tboot = databoot(:,1)+databoot(:,2);
        fwd_func = @(x)forward_TD(x,dynx,Tboot,Tboot(1),ctHb,Fs);
        xfitboot(i,:) = fmincon(@(x0)costfunc(fwd_func,x0,databoot,data_err),x0,[],[],[],[],lb,ub,[],options); 
    end
    fprintf('\n');
    xfit_s = std(xfitboot);
end

%% Cost function
function chi2r = costfunc(fwd_func,x,data,data_err) 
    [O,D] = feval(fwd_func,x);
    chi2_O = sum(((O-data(:,1))./data_err(:,1)).^2);
    chi2_D = sum(((D-data(:,2))./data_err(:,2)).^2);
    chi2r = (chi2_O+chi2_D)./(2*length(O)-length(x)-1);   
end