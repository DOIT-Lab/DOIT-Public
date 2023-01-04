function [outp] = DPF_DSF_calc(rs, rd, muaLamb, muspLamb, calcmode,...
    fmod, nin)
% [outp] = DPF_DSF_calc(rs, rd, muaLamb, muspLamb, calcmode, fmod, nin)
% Updated by Giles Blaney Spring 2022
% Inputs:
%   rs          - Pencil beam source positions in [x, y, z] with direction 
%                   vector of [0, 0, 1].
%                   Size: [numSD, 3]; Units: mm
%   rd          - Detector positions in [x, y, z].
%                   Size: [numSD, 3]; Units: mm
%   muaLamb     - Vector of absorption coefficients.
%                   Size: [numOptProp, 1] OR [1, numOptProp]; Units: 1/mm
%   muspLamb    - Vector of reduced scattering coefficients.
%                   Size: [numOptProp, 1] OR [1, numOptProp]; Units: 1/mm
%   calcmode    - String with four possible values:
%                   - 'DPF_I': To calculate the intensity DPF
%                   - 'DPF_Ph' or 'DPF_P': To calculate the phase DPF
%                   - 'DSF_I': To calculate the intensity DSF
%                   - 'DSF_Ph' or 'DSF_P': To calculate the phase DSF
%   fmod        - OPTIONAL, default=140.625e6 Hz
%                   Modulation frequency
%                   Size: scalar; Units: Hz
%   nin         - OPTIONAL, default=1.4
%                   Index of refraction inside the medium (outside assumed
%                   to be 1).
%                   Size: scalar; Units: unitless
% Outputs:
%   outp        - Value requested by calcmode input.
%                   Units: unitless
%                   - For DSF or DPF if numSD=1
%                       Size: same as muaLamb input
%                   - For DPF if numSD>1
%                       Size: [numOptProp, numSD]

    %% Setup
    if nargin<6
        fmod=1.40625e8; %Hz
        nin=1.4;
    end
    nout=1;
    
    rho=(sqrt(sum((rs-rd).^2, 2))).';
    rhoAvg=mean(rho);

    if strcmp(calcmode(1:3), 'DPF') && length(rho)~=1
        outp=NaN(length(muaLamb), length(rho));
    else
        outp=NaN(size(muaLamb));
    end
    
    if calcmode(end)=='C'
        fmod=0;
        calcmode(end)='I';
    end
    
    optProp.nin=nin;
    optProp.nout=nout;

    %% Calc DPF or DSF
    for lamInd=1:length(muaLamb)
        optProp.mua=muaLamb(lamInd);
        optProp.musp=muspLamb(lamInd);
        
        Lac=NaN(size(rho));
        Lph=NaN(size(rho));
        rs_iso=rs+[0, 0, 1./optProp.musp];
        for i=1:length(rho)
            L=complexTotPathLen(rs_iso(i, :), rd(i, :),...
                2*pi*fmod, optProp);
            Lac(i)=real(L);
            Lph(i)=imag(L);
        end        
        
        switch calcmode
            case 'DPF_I'
                if strcmp(calcmode(1:3), 'DPF') && length(rho)~=1
                    outp(lamInd, :)=Lac./rho;
                else
                    outp(lamInd)=Lac./rho;
                end

            case {'DPF_Ph', 'DPF_P'}
                if strcmp(calcmode(1:3), 'DPF') && length(rho)~=1
                    outp(lamInd, :)=Lph./rho;
                else
                    outp(lamInd)=Lph./rho;
                end

            case 'DSF_I'
                outp(lamInd)=sum((rho-rhoAvg).*Lac, 2)/...
                    (length(rho)*var(rho, 1));

            case {'DSF_Ph', 'DSF_P'}
                outp(lamInd)=sum((rho-rhoAvg).*Lph, 2)/...
                    (length(rho)*var(rho, 1));

            otherwise
                error('Unknown calcmode');
        end
    end
end