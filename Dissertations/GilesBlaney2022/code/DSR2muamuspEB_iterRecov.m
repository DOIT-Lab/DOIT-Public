function [mua, musp, iter]=DSR2muamuspEB_iterRecov(rhos, RR, opts)
% [mua, musp, iter]=DSR2muamuspEB_iterRecov(rhos, RR, opts)
% Giles Blaney Fall 2020
% Assumes extrapolated-boundary condition
% Expects DS set made of 4 SD measurements
%   Inputs:
%       rhos     - 1 x 4 array of source-detector distances in the
%                  following format: [S1, L1, S2, L2] (mm)            
%       RR       - Complex relectance data in a 1 x 4 array containing data
%                  in the following format: [S1, L1, S2, L2]
%       opts     - (Optional) Structure containing options as feilds:
%                  - mua0: mua to use at for first iteration (1/mm)
%                  - musp0: mua to use at for first iteration (1/mm)
%                  - mueff_tol: Change in mueff stoping criteria (1/mm)
%                    abs(mueff(n)-mueff(n-1))<mueff_tol => Stop
%                  - n_max: Maximum iteration number stoping criteria
%                    n>=n_max => Stop
%                  - omega: Modulation frequency (rad/sec) 
%                  - v: Speed of light in tissue (mm/sec)
%                  - nin: Index of refraction inside medium
%                  - nout: Index of refraction outside medium
%   Output:
%       mua      - mua found (1/mm)
%       mua      - musp found (1/mm)
%       iter     - Structure containing information about iterations as
%                  feilds:
%                  - n: Number of iterations
%                  - mueff_all: History of mueff during iterations
%                  - mua_all: History of mua during iterations
%                  - musp_all: History of musp during iterations

    %% Parse Input
    if nargin<=2
        mu0Bool=true;
        mueff_tol=1e-4; %1/mm
        n_max=10;
        omega=2*pi*140.625e6; %rad/sec
        nin=1.4;
        v=2.99792458e11/nin; %mm/sec
        nout=1;
    else
        if isfield(opts, 'mua0')
            mua0=opts.mua0;
            musp0=opts.musp0;
            mu0Bool=false;
        else
            mu0Bool=true;
        end
        mueff_tol=opts.mueff_tol;
        n_max=opts.n_max;
        omega=opts.omega;
        v=opts.v;
        nin=opts.nin;
        nout=opts.nout;
    end
    
    if length(rhos)<4
        rhos=[rhos, rhos];
    end
    
    rho1=rhos(1:2);
    rho2=rhos(3:4);
    R1=RR(1:2);
    R2=RR(3:4);
    
    if angle(R1(1))>angle(R1(2))
        R1(2)=abs(R1(2))*exp(1i*(angle(R1(2))+2*pi));        
    end
        
    if angle(R2(1))>angle(R2(2))
        R2(2)=abs(R2(2))*exp(1i*(angle(R2(2))+2*pi));        
    end
    
    %% Find mua0
    if mu0Bool
        I1=abs(R1);
        I2=abs(R2);
        P1=angle(R1);
        P2=angle(R2);
        
        SSI10=(log(rho1(2)^2*I1(2))-log(rho1(1)^2*I1(1)))/...
            (rho1(2)-rho1(1));
        SSI20=(log(rho2(2)^2*I2(2))-log(rho2(1)^2*I2(1)))/...
            (rho2(2)-rho2(1));
        
        SSP10=wrapToPi(P1(2)-P1(1))/...
            (rho1(2)-rho1(1));
        SSP20=wrapToPi(P2(2)-P2(1))/...
            (rho2(2)-rho2(1));

        DSI0=(SSI10+SSI20)/2;
        DSP0=(SSP10+SSP20)/2;

        mua0=(omega/(2*v))*(DSP0./DSI0-DSI0./DSP0);
        musp0=(DSI0.^2-DSP0.^2)./...
            (3*mua0)-mua0;
    end
    
    %% Iterative mua Recovery
    n=1;
    stopBool=false;
    mua=mua0; %1/mm
    musp=musp0; %1/mm
    mueff0=sqrt(3*(mua0-1i*omega/v)*(musp0+mua0));
    A=n2A(nin, nout);
    mueff=mueff0;
    while ~stopBool
        z0=1/musp(n);
        zb=-2*A/(3*(musp(n)+mua(n)));
        z0p=-z0+2*zb;
        r1_1=sqrt(rho1.^2+z0.^2);
        r1_2=sqrt(rho2.^2+z0.^2);
        r2_1=sqrt(rho1.^2+z0p.^2);
        r2_2=sqrt(rho2.^2+z0p.^2);
        
        C1_1=z0.*(1./r1_1+mueff(n))./(r1_1.^2);
        C1_2=z0.*(1./r1_2+mueff(n))./(r1_2.^2);
        C2_1=-z0p.*(1./r2_1+mueff(n))./(r2_1.^2);
        C2_2=-z0p.*(1./r2_2+mueff(n))./(r2_2.^2);
        
        y_1=log(4*pi*R1./...
            (C1_1+C2_1.*exp(mueff(n).*(r1_1-r2_1))));
        y_2=log(4*pi*R2./...
            (C1_2+C2_2.*exp(mueff(n).*(r1_2-r2_2))));

        SSR1=diff(y_1)/diff(r1_1);
        SSR2=diff(y_2)/diff(r1_2);
        
        DSR=(SSR1+SSR2)/2;

        n=n+1;

        mueff(n)=-DSR;
        
        a=real(mueff(n));
        b=imag(mueff(n));
        mutp=-2*a*b*v/(3*omega);
        mua(n)=(a^2-b^2)/(3*mutp);
        musp(n)=(mutp-mua(n));
        
        if n>=n_max
            warning('Max iterations reached.');
            stopBool=true;
        elseif abs(mueff(end)-mueff(end-1))<mueff_tol
            stopBool=true;
        end
    end
    
    %% Package Output
    iter.n=n;
    iter.mueff_all=mueff;
    iter.mua=mua;
    iter.musp=musp;
    
    mua=mua(end);
    musp=musp(end);
    
end