function [O,D] = forward_TD(x,dynx,T,T0,ctHb,Fs)
    % Time-domain coherent hemodynamics spectroscopy model
    % Inputs:
    %       x  -   CHS baseline parameters
    %              x = [tc tv F(c)CBV0c/CBV0 fc k]. Units are [s s - Hz -].
    %       dynx - CHS dynamic parameters: blood volume contributions from
    %               the arteries, capillaries, and veins to the total blood volume
    %              dynx = [CBVa/CBV, CBVc/CBV, CBVv/CBV].
    %       T   -   Temporal trace of absolute totalhemoglobin concentration. 
    %               Unit is microM.
    %       T0  -   Absolute totalhemoglobin concentration. Unit is microM.
    %       ctHb -  hemoglobin concentration in blood. Unit is microM/unit of blood volume
    %       Fs  -   Sampling rate. Unit is Hz
    % Outputs:
    %       O   -   Temporal trace of absolute oxyhemoglobin concentration. 
    %               Unit is microM.
    %       D   -   Temporal trace of absolute deoxyhemoglobin concentration. 
    %               Unit is microM.

    % Written By Jana Kainerstorfer
    % Modified by Kristen Tgavalekos
    % Last modified by Thao Pham (February 2020)
    
    L = length(T);
    t  = 0:1/Fs:2000;
    t  = t-max(t)/2;
    tstart = find(t>=0,1);
    
    % Parameters of the CHS model
    tc  = x(1); % capillary blood transit time (s)
    tv  = x(2); % venous blood transit time (s)
    Fpc = x(3); % capillary baseline blood volume ratio (F(c)CBV0(c)/CBV0)
    fc  = x(4); % cut-off frequency for autoregulation high-pass filter (Hz)
    k   = x(5); % inverse of the modified Grubb exponent

    a = dynx(1); % CBVa/CBV
    c = dynx(2); % CBVc/CBV
    v = dynx(3); % CBVv/CBV
    
    % Fixed paramters
    pa = 0.3; % arterial baseline blood volume ratio (CBV0(a)/CBV0)
    alpha = 0.8; % rate constant of oxygen diffusion (1/s)
    Sa  = 0.98; % arterial blood saturation
    Sc  = Sa*(1/tc/alpha)*(1-exp(-alpha*tc)); % capillary blood saturation
    Sv  = Sa*exp(-alpha*tc); % venous blood saturation

    % Static blood volume
    pv     = 1-(pa+Fpc); % venous baseline blood volume ratio (CBV0(v)/CBV0)
    CBV0   = T0/ctHb; % baseline blood volume 
    CBV0a  = pa*CBV0; % arterial baseline blood volume (CBV0(a))
    FCBV0c = Fpc*CBV0; % capillary baseline blood volume (F(c)CBV0(c))
    CBV0v  = pv*CBV0; % venous baseline blood volume (CBV0(v))

    % Frequency response
    h_RC_LP = (exp(1)/tc)*exp(-t*exp(1)/tc); % RC lowpass
    h_RC_LP(t<0) = 0;
    h_G_LP = (1/0.6/(tc+tv))*... % Gaussian lowpass
        exp((-pi*(t-0.5*(tc+tv)).^2)/((0.6*(tc+tv))^2));
    tau =  1/(fc*2*pi);
    h_RC_HP = -(1/tau)*exp(-t/tau); % RC highpass
    h_RC_HP(t<0) = 0;
    
    % cerebral blood volume change (cbv)   
    cbv = zeros(1,size(t,2));
    cbv(tstart:L+tstart-1) = (T-T0)/T0;
    CBV = CBV0*cbv;    
    CBVa = a*CBV; % arterial blood volume change
    CBVc = c*CBV; % capillary blood volume change
    CBVv = v*CBV; % venous blood volume change

    % cerebral blood flow change (cbf)
    cbftmp = (1/Fs)*conv(h_RC_HP,cbv);
    cbftmp = cbftmp(tstart:end);
    cbftmp = cbftmp(1:size(t,2));
    cbf    = k*(cbv+cbftmp);

    %% Time-dependent expressions for absolute hemoglobin concentrations
    % O 
    O0  = ctHb*(Sa*CBV0a+Sc*FCBV0c+Sv*CBV0v); % absolute baseline oxyHb (microM)
    dOv = ctHb*(Sa*CBVa+Sc*CBVc+Sv*CBVv); % oxyHb response due to blood volume change (microM)
    dOf = ctHb*(1/Fs)*...% oxyHb response due to blood flow change (microM)
        conv((Sc/Sv)*(Sc-Sv)*FCBV0c*h_RC_LP+(Sa-Sv)*CBV0v*h_G_LP,cbf);
    dOf = dOf(tstart:end);
    dOf = dOf(1:size(t,2));
    O   = O0+dOv+dOf;
    O   = O(tstart:tstart+L-1)';

    % D 
    D0  = ctHb*((1-Sa)*CBV0a+(1-Sc)*FCBV0c+(1-Sv)*CBV0v); % absolute baseline deoxyHb (microM)
    dDv = ctHb*((1-Sa)*CBVa+(1-Sc)*CBVc+(1-Sv)*CBVv);% deoxyHb response due to blood volume (microM)
    dDf = -ctHb*(1/Fs)*...% deoxyHb response due to blood flow (microM)
        conv((Sc/Sv)*(Sc-Sv)*FCBV0c*h_RC_LP+(Sa-Sv)*CBV0v*h_G_LP,cbf);
    dDf = dDf(tstart:end);
    dDf = dDf(1:size(t,2));
    D   = D0+dDv+dDf;
    D   = D(tstart:tstart+L-1)';
end

