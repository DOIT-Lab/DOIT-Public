function [cohThresh, t, f, n, eltime, alpha, coi]=...
    genWAVthresh(rtime, tTot, fs, alpha, svBool)
% [cohThresh, t, f, n, eltime]=genWAVthresh(rtime, tTot, fs, alpha, svBool)
% By Giles Blaney Fall 2018
%   Inputs:
%       rtime     - Approximate time to run in seconds
%       tTot      - Total length of protocol in seconds
%       fs        - Sampling frequency in Hz
%       alpha     - Significance level
%       svBool    - Boolean to save threshold map
%   Outputs:
%       cohThresh - Wavelet coherence threshold map
%       t         - Time vector seconds
%       f         - Frequency vector in Hz
%       n         - Number of samples run
%       eltime    - Time to run in seconds
%       alpha     - Significance level

    % Set Defaults
    if nargin<=0
        rtime=60; %sec
    end
    if nargin<=1
        tTot=10*60; %sec
    end
    if nargin<=2
        fs=9.9306; %Hz
    end
    if nargin<=3
        alpha=0.05;
    end
    if nargin<=4
        svBool=true;
    end

    nsamp=floor(tTot*fs); % Number of samples in time
    t=(0:(1/fs):((nsamp-1)/fs))'; % Create time vector (sec)
    
    % Run coherence on signals of correct length to get frequency vector,
    % cone of influence vector, and example coherence array
    sig1=ones(nsamp, 1);
    sig2=ones(nsamp, 1);
    [Xtemp, ~, ~, ~, f, coi]=wcoherence(sig1, sig2, fs);
    
    % Set size of A set and B set
    Nkp=100; % Number of samples to keep each iteration
    AN=floor(Nkp*0.5); % A set is half of all kept samples
    BN=ceil(Nkp*0.5); % B set is other half of all kept samples
    Ainit=round(Nkp*(1-alpha)); % Size of A in initialization step
    Binit=Nkp-Ainit; % Size of B in initialization step

    % Get expected size of coherence array and initialize AB array and
    % threshold array of NaN
    N1=size(Xtemp, 1);
    N2=size(Xtemp, 2);
    AB=NaN*ones(N1, N2, Nkp+1);
    cohThresh=NaN*ones(N1, N2);
    
    % Clear variables that will not be used later
    clear tTot Xtemp N1 N2;
    
    tic; % Start timer
    n=0; % Initialize iteration number
    eltime=toc; % Initialize current time
    % Run while current time is less than run time or iteration number is
    % less then the amount of samples needed to initialize AB
    while eltime<rtime || n<Nkp
        n=n+1; % Increment iteration number
        % Generate Gaussian random signals
        sig1=randn(nsamp, 1);
        sig2=randn(nsamp, 1);
        
        % Calculate coherenace of random signals
        X_now=wcoherence(sig1, sig2, fs);
        
        if n<=Nkp % If AB is not initialized
            AB(:, :, n)=X_now; % Populate AB with current sample

            %% Initialization of AB
            if n==Nkp % If at last step of AB initialization
                % Set number of samples in Aall set to the size of A at
                % initialization
                An=ones(size(X_now))*Ainit;

                % Set all values too small to be in A to NaN and sort
                AB=sort(AB, 3);
                AB(:, :, 1:(BN-Binit))=NaN;
                AB=sort(AB, 3);

                % Initialize threshold array with the first page of B
                cohThresh=AB(:, :, AN+1);
            end
        else % If AB is initialized
            %% Add New Sample
            % Place the current sample in the last page of AB
            AB(:, :, end)=X_now;

            % Find the number of new samples that belong in A, increment
            % the number of samples in Aall if this is the case, and sort
            % AB with the new sample
            Aadd=AB(:, :, end)<=cohThresh;
            An=An+double(Aadd);
            AB=sort(AB, 3);

            % If the current sample was placed in A set the first page of A
            % to NaN and sort
            Amin=AB(:, :, 1);
            Amin(Aadd)=NaN;
            AB(:, :, 1)=Amin;
            AB=sort(AB, 3);

            %% Shift Samples
            % Find if the number of samples in Aall is too large or too
            % small given the current divide between A and B
            Asmall=An<((1-alpha)*n);
            Abig=An>((1-alpha)*n);

            % If Aall is too small set the first page of A to NaN and
            % increment Aall, if Aall is too big set the last page of B to
            % -Inf, then sort AB
            Amin=AB(:, :, 1);
            Bmax=AB(:, :, end);
            Amin(Asmall)=NaN;
            An(Asmall)=An(Asmall)+1;
            Bmax(Abig)=-Inf;
            An(Abig)=An(Abig)-1;
            AB(:, :, 1)=Amin;
            AB(:, :, end)=Bmax;   
            AB=sort(AB, 3);

            % Set the threshold to the first element in B
            cohThresh=AB(:, :, AN+1);

            %% Check for Problem Elements
            % Check if NaN or -Inf where placed in the threshold estimate,
            % if so replace with the max sample value
            thresh_bad=or(isinf(cohThresh), isnan(cohThresh));
            ABmax=max(AB, [], 3, 'omitnan');
            cohThresh(thresh_bad)=ABmax(thresh_bad);
        end
        
        eltime=toc; % Get current time
        
        % If iteration number is multiple of 100 print status
        if mod(n, 100)==0
            fprintf('Sample %d\t%.1f sec / %.0f sec\n', n, eltime, rtime);
        end
    end
    fprintf('Loop Done\n');
    
    if svBool % If desired save output
        % Create string of sampling frequecy replacing . with o
        fsString=num2str(fs);
        fsString=replace(fsString, '.', 'o');
        
        % Find file path for this function
        funPathStr=which('genWAVthresh.m');
        funDir=dir(funPathStr);
        
        % Construct filename to be saved
        svString=[funDir.folder '\WAVcohThresh_' fsString '_' ...
            num2str(nsamp) '_' num2str(n) '.mat'];
        
        % Save
        save(svString, 'cohThresh', 't', 'f', 'n', 'eltime', 'alpha', 'coi');
        
        fprintf('Saving Done\n');
    end
    
end