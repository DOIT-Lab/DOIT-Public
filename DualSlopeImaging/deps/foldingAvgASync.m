function [tFold, xFold, sFold, XX] = foldingAvgASync(t, x, t0, dtn, dtp,...
    opts)
% [tFold, xFold, sFold, XX] = foldingAvgASync(t, x, t0, dtn, dtp, dtBool)
% Giles Blaney Updated Summer 2021
% Inputs:
%           -    t: Time vector
%           -    x: Data vector
%           -   t0: Vector of event times
%           -  dtn: Positive time interval after t0 to fold
%           -  dtp: Negative time interval before t0 to fold
%           - opts: Struct of options
%               - DTbools(1): Boolean to apply linear detrending for each
%                             event
%               - DTbools(2): Boolean to apply linear detrending over all
%                             events
%               -        nDT: Number of samples at begining and end of fold
%                             interval for event detrending 
%                             (for DTbools(1))
% Outputs:
%           - tFold: Folded time vector
%           - xFold: Folded data vector
%           - sFold: Vector of standard deviations of folded data
%           -    XX: Time x event matrix of data pre-average
    
    if nargin<=5
        opts.DTbools=[false, false];
        opts.nDT=10;
    end
    if ~isstruct(opts) % For legacy input (6th input DTbools)
        opts.DTbools=opts;
    end
    
    if size(t, 2)~=1
        t=t.';
    end
    if length(t)~=size(x, 1)
        x=x.';
    end
    
    XX=[];
    for i=1:length(t0)
        tmin=t0(i)-dtn;
        tmax=t0(i)+dtp;
        
        xTemp=x(and(t>tmin, t<=tmax)); % changed
        
        if opts.DTbools(1)
            xTemp=xTemp-mean(xTemp(1:opts.nDT), 'omitnan');
            
            xTemp=xTemp-mean(xTemp((end-opts.nDT):end), 'omitnan').*...
                (0:(length(xTemp)-1)).'/(length(xTemp)-1);
        end
        
        if i==1
            XX=xTemp;
        else
            while length(xTemp)>size(XX, 1)
                XX=[XX; NaN(1, size(XX, 1))];
            end
            while length(xTemp)<size(XX, 1)
                xTemp=[xTemp; NaN];
            end
            XX=[XX, xTemp];
        end
    end
    
    xFold=mean(XX, 2, 'omitnan');
    sFold=std(XX, [], 2, 'omitnan');
    tFold=timeAxis(1/median(diff(t)), length(xFold));
    
    if length(tFold)>opts.nDT
        xFold=xFold-mean(xFold(1:opts.nDT), 'omitnan');
    else
        xFold=xFold-mean(xFold, 'omitnan');
    end
    
    tFold=tFold-dtn;
end

