function [T] = Tcuv(rSrc, rDet, L, optProp, opts)
    
    % Src is in +y dir 
    % rSrc & rDet = [x, z]
    % L = [x, y, z]

    %% Setup
    if nargin<=4
        omega=2*pi*140.625e6; %rad/sec
        BC='EBC';

        lMax=3;
        mMax=3;
        nMax=1;
        
        warning(['Default opts used']);
    else
        omega=opts.omega;
        BC=opts.BC;
        
        lMax=opts.lMax;
        mMax=opts.mMax;
        nMax=opts.nMax;
    end
    
    mua=optProp.mua;
    musp=optProp.musp;
    nin=optProp.nin;
    nout=optProp.nout;

    c=2.99792458e11; %mm/sec
    v=c/nin;
    
    switch BC
        case 'EBC'
            A=n2A(nin, nout);
        case 'ZBC'
            A=0;
        otherwise
            error('Unknown BC');
    end
    
    %Source position
    xu=rSrc(1);
    y0=1/musp;
    zu=rSrc(2);

    D=1/(3*musp);
    ue0=sqrt((mua*v-1i*omega)/(v*D));% with this definition the phase
    % increases as the s-d separation increases

    yb=2*A*D; %extrapolated length

    %Detector position
    xd=rDet(1);
    yd=L(2);
    zd=rDet(2);
    
    %% Calc
    l=-lMax:lMax;
    m=-mMax:mMax;
    n=-nMax:nMax;

    % Two ways of defining the source positions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) This way is less elegant but we can have different number of 
    % sources in the three directions. This is useful in the cuvette 
    % geometry where the number of sources along the longest side affect 
    % little light propagation. There are
    % 8*(length(l) x length(m) x length(n)) different point sources
    % (including positive and negative)
    
    x1l=2*L(1)*l+4*l*yb+xu; %positive point sources
    x2l=2*L(1)*l+(4*l-2)*yb-xu; %negative point sources
    
    y1m=2*L(2)*m+4*m*yb+y0; %positive point sources
    y2m=2*L(2)*m+(4*m-2)*yb-y0; %negative point sources
    
    z1n=2*L(3)*n+4*n*yb+zu; %positive point sources
    z2n=2*L(3)*n+(4*n-2)*yb-zu; %negative point sources

    % ---------------------------------------------------------------------
    % 2) This way is more elegant but we require the same number of sources
    % in the three directions.
    % 
    % Definition of the point source positions: there are
    % 8*(length(l))^(3) different point sources (including positive and
    % negative)
    % 
    % N=[l' m' n'];
    % P1=2*N.*L+4*N*yb+[xu, y0, zu]; %positive point sources (length(l)x3))
    % P2=2*N.*L+(4*N-2)*yb-[xu, y0, zu]; %negative point sources
    % 
    % Coordinates of positive sources
    % x1l=P1(:,1);
    % y1m=P1(:,2);
    % z1n=P1(:,3);
    % 
    % x2l=P2(:,1);
    % y2m=P2(:,2);
    % z2n=P2(:,3);
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Xd,Yd,Zd] = meshgrid(...
        xd*ones(length(l),1), yd*ones(length(m),1), zd*ones(length(n),1));
    
    [X1,Y1,Z1] = meshgrid(x1l,y1m,z1n);
    r1=sqrt((X1-Xd).^2+(Y1-Yd).^2+(Z1-Zd).^2);
    
    [X1,Y2,Z1] = meshgrid(x1l,y2m,z1n);
    r2=sqrt((X1-Xd).^2+(Y2-Yd).^2+(Z1-Zd).^2);
    
    [X1,Y1,Z2] = meshgrid(x1l,y1m,z2n);
    r3=sqrt((X1-Xd).^2+(Y1-Yd).^2+(Z2-Zd).^2);
    
    [X1,Y2,Z2] = meshgrid(x1l,y2m,z2n);
    r4=sqrt((X1-Xd).^2+(Y2-Yd).^2+(Z2-Zd).^2);
    
    [X2,Y1,Z1] = meshgrid(x2l,y1m,z1n);
    r5=sqrt((X2-Xd).^2+(Y1-Yd).^2+(Z1-Zd).^2);
    
    [X2,Y2,Z1] = meshgrid(x2l,y2m,z1n);
    r6=sqrt((X2-Xd).^2+(Y2-Yd).^2+(Z1-Zd).^2);
    
    [X2,Y1,Z2] = meshgrid(x2l,y1m,z2n);
    r7=sqrt((X2-Xd).^2+(Y1-Yd).^2+(Z2-Zd).^2);
    
    [X2,Y2,Z2] = meshgrid(x2l,y2m,z2n);
    r8=sqrt((X2-Xd).^2+(Y2-Yd).^2+(Z2-Zd).^2);
    
    T_mat=1/(4*pi)*((L(2)-Y1).*(ue0+1./r1).*exp(-ue0*r1)./r1.^2-...
        (L(2)-Y2).*(ue0+1./r2).*exp(-ue0*r2)./r2.^2-...
        (L(2)-Y1).*(ue0+1./r3).*exp(-ue0*r3)./r3.^2+...
        (L(2)-Y2).*(ue0+1./r4).*exp(-ue0*r4)./r4.^2-...
        (L(2)-Y1).*(ue0+1./r5).*exp(-ue0*r5)./r5.^2+...
        (L(2)-Y2).*(ue0+1./r6).*exp(-ue0*r6)./r6.^2+...
        (L(2)-Y1).*(ue0+1./r7).*exp(-ue0*r7)./r7.^2-...
        (L(2)-Y2).*(ue0+1./r8).*exp(-ue0*r8)./r8.^2);
    
    T_mat=reshape(T_mat,length(m),length(l),length(n));
    T=squeeze(sum(T_mat,[1 2 3]));
end