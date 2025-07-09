function [mua_RelErr, musp_RelErr, ...
    dmuaDSI_RelErr, dmuaDSP_RelErr, ...
    optProp_rec] = DSpos2err(...
    rsPen_act, rd_act, rsPen_ass, rd_ass, optProp_act, dmua_act, omega, ...
    recMeth)
arguments
    rsPen_act (2,3) double;
    rd_act (2,3) double;

    rsPen_ass (2,3) double;
    rd_ass (2,3) double;
    
    optProp_act struct;
    % Contains: mua, musp, nin, nout
    
    dmua_act (1,1) double = 0.0001; %1/mm
    omega (1,1) double = 140.625e6*2*pi; %rad/s

    recMeth string = "iter"; %"iter" or "lin"
end
    
    %% Compute Forward R
    rs_act=rsPen_act+[0, 0, 1/optProp_act.musp];
    
    R1A=complexReflectance(rs_act(1, :), rd_act(1, :), ...
        omega, optProp_act);
    R1B=complexReflectance(rs_act(1, :), rd_act(2, :), ...
        omega, optProp_act);
    R2A=complexReflectance(rs_act(2, :), rd_act(1, :), ...
        omega, optProp_act);
    R2B=complexReflectance(rs_act(2, :), rd_act(2, :), ...
        omega, optProp_act);
    RR=[R1A, R1B, R2B, R2A];

    %% Recover Absolute mua and musp
    rho1A_ass=vecnorm(rsPen_ass(1, :)-rd_ass(1, :));
    rho1B_ass=vecnorm(rsPen_ass(1, :)-rd_ass(2, :));
    rho2A_ass=vecnorm(rsPen_ass(2, :)-rd_ass(1, :));
    rho2B_ass=vecnorm(rsPen_ass(2, :)-rd_ass(2, :));
    rhos_ass=[rho1A_ass, rho1B_ass, rho2B_ass, rho2A_ass];
    
    optsRec.omega=omega;
    optsRec.nin=optProp_act.nin;
    optsRec.nout=optProp_act.nout;
    [mua_tmp, musp_tmp, iter_tmp]= ...
        DSR2muamuspEB_iterRecov(rhos_ass, RR, optsRec);

    switch recMeth
        case "iter"
            mua_rec=mua_tmp;
            musp_rec=musp_tmp;
        case "lin"
            mua_rec=iter_tmp.mua(1);
            musp_rec=iter_tmp.musp(1);
        otherwise
            error('Unknown optical prop rec method');
    end
    
    optProp_rec.mua=mua_rec;
    optProp_rec.musp=musp_rec;
    optProp_rec.nin=optsRec.nin;
    optProp_rec.nout=optsRec.nout;
    
    mua_RelErr=(optProp_rec.mua-optProp_act.mua)/optProp_act.mua;
    musp_RelErr=(optProp_rec.musp-optProp_act.musp)/optProp_act.musp;

    %% Compute <L> from Recovered mua and musp, and assumed rhos
    L1A= complexTotPathLen([0, 0, 1/optProp_rec.musp], [rho1A_ass, 0, 0], ...
        omega, optProp_rec);
    L1B= complexTotPathLen([0, 0, 1/optProp_rec.musp], [rho1B_ass, 0, 0], ...
        omega, optProp_rec);
    L2A= complexTotPathLen([0, 0, 1/optProp_rec.musp], [rho2A_ass, 0, 0], ...
        omega, optProp_rec);
    L2B= complexTotPathLen([0, 0, 1/optProp_rec.musp], [rho2B_ass, 0, 0], ...
        omega, optProp_rec);

    %% Compute Forward dR from small dmua
    optProp1_act=optProp_act;
    optProp1_act.mua=optProp_act.mua+dmua_act;

    R1A_1=complexReflectance(rs_act(1, :), rd_act(1, :), ...
        omega, optProp1_act);
    R1B_1=complexReflectance(rs_act(1, :), rd_act(2, :), ...
        omega, optProp1_act);
    R2A_1=complexReflectance(rs_act(2, :), rd_act(1, :), ...
        omega, optProp1_act);
    R2B_1=complexReflectance(rs_act(2, :), rd_act(2, :), ...
        omega, optProp1_act);

    %% Recover dmua from all other recovered and assumed values
    %DSI
    SSI1_0=(log(rho1B_ass^2*abs(R1B))-log(rho1A_ass^2*abs(R1A)))/...
        (rho1B_ass-rho1A_ass);
    SSI2_0=(log(rho2A_ass^2*abs(R2A))-log(rho2B_ass^2*abs(R2B)))/...
        (rho2A_ass-rho2B_ass);
    SSI1_1=(log(rho1B_ass^2*abs(R1B_1))-log(rho1A_ass^2*abs(R1A_1)))/...
        (rho1B_ass-rho1A_ass);
    SSI2_1=(log(rho2A_ass^2*abs(R2A_1))-log(rho2B_ass^2*abs(R2B_1)))/...
        (rho2A_ass-rho2B_ass);

    DSFI1=(real(L1B)-real(L1A))/(rho1B_ass-rho1A_ass);
    DSFI2=(real(L2A)-real(L2B))/(rho2A_ass-rho2B_ass);
    
    dmuaDSI=-((SSI1_1-SSI1_0)+(SSI2_1-SSI2_0))/...
        (DSFI1+DSFI2);

    dmuaDSI_RelErr=(dmuaDSI-dmua_act)/dmua_act;

    %DSP
    SSP1_0=(angle(R1B)-angle(R1A))/...
        (rho1B_ass-rho1A_ass);
    SSP2_0=(angle(R2A)-angle(R2B))/...
        (rho2A_ass-rho2B_ass);
    SSP1_1=(angle(R1B_1)-angle(R1A_1))/...
        (rho1B_ass-rho1A_ass);
    SSP2_1=(angle(R2A_1)-angle(R2B_1))/...
        (rho2A_ass-rho2B_ass);

    DSFP1=(imag(L1B)-imag(L1A))/(rho1B_ass-rho1A_ass);
    DSFP2=(imag(L2A)-imag(L2B))/(rho2A_ass-rho2B_ass);
    
    dmuaDSP=-((SSP1_1-SSP1_0)+(SSP2_1-SSP2_0))/...
        (DSFP1+DSFP2);

    dmuaDSP_RelErr=(dmuaDSP-dmua_act)/dmua_act;
    
end