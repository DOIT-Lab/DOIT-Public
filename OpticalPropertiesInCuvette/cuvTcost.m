function [cost] = cuvTcost(mu, Tmeas, Terr, rSrc, rDet, opts)
    
    %% DT
    nDet=size(rDet, 1);
    if nDet~=2
        error('Only use two detectors');
    end
    T=NaN(nDet, 1);
    optProp.mua=mu(1);
    optProp.musp=mu(2);
    optProp.nin=opts.nin;
    optProp.nout=opts.nout;
    L=opts.L;
    for i=1:nDet
        T(i)=Tcuv(rSrc, rDet(i, :), L, optProp, opts);
    end
    
    logAmpDiff=log(abs(T(2)))-log(abs(T(1)));
    phDiff=angle(T(2))-angle(T(1));
    
    if size(Tmeas, 2)==1
        logAmpDiff_meas=log(abs(Tmeas(2)))-log(abs(Tmeas(1)));
        phDiff_meas=angle(Tmeas(2))-angle(Tmeas(1));
    else
        logAmpDiff_meas1=log(abs(Tmeas(2, 1)))-log(abs(Tmeas(1, 1)));
        phDiff_meas1=angle(Tmeas(2, 1))-angle(Tmeas(1, 1));

        logAmpDiff_meas2=log(abs(Tmeas(2, 2)))-log(abs(Tmeas(1, 2)));
        phDiff_meas2=angle(Tmeas(2, 2))-angle(Tmeas(1, 2));

        logAmpDiff_meas=(logAmpDiff_meas1+logAmpDiff_meas2)/2;
        phDiff_meas=(phDiff_meas1+phDiff_meas2)/2;
    end
    
    logAmpDiff_err=(sqrt((abs(Terr(1))/abs(Tmeas(1)))^2+...
        (abs(Terr(2))/abs(Tmeas(2)))^2))/sqrt(2);
    phDiff_err=(sqrt(angle(Terr(1))^2+angle(Terr(2))^2))/sqrt(2);
    
    cost=opts.kappa*((logAmpDiff_meas-logAmpDiff)/logAmpDiff_err)^2+...
        ((phDiff_meas-phDiff)/phDiff_err)^2;
end