function cbf_cmro2 = cbf_cmro2_dyn(x,dynx,T0,dO,dD,Fs)
    % coherent hemodynamics spectroscopy model to calculate cbf-cmro2
    % Inputs:
    %       x  -   CHS baseline parameters
    %              x = [tc tv F(c)CBV0c/CBV0]. Units are [s s -].
    %       dynx - CHS dynamic parameters: blood volume contributions from
    %               the arteries, capillaries, and veins to the total blood volume
    %              dynx = [CBVa/CBV, CBVc/CBV, CBVv/CBV].
    %       T0  -   Absolute totalhemoglobin concentration. Unit is microM.
    %       dO  -   Temporal trace of relative oxyhemoglobin concentration. 
    %               Unit is microM.
    %       dD  -   Temporal trace of relative deoxyhemoglobin concentration. 
    %               Unit is microM.
    %       Fs  -   Sampling rate. Unit is Hz
    % Outputs:
    % cbf-cmro2   -   Temporal trace of cbf-cmro2

    % Written By Jana Kainerstorfer
    % Modified by Kristen Tgavalekos
    % Last modified by Thao Pham (February 2020)
    
    % Setups
    dT_t = dO + dD;
    cbv_t = dT_t/T0;

    L = length(dO);
    n = 2^nextpow2(L);
    f = (0:n-1)*(Fs/n);
    dO_w = (1/Fs)*fft(dO,n);
    dD_w = (1/Fs)*fft(dD,n);
    cbv_w = (1/Fs)*fft(cbv_t,n);

    ctHb = 2300;
    CBV0 = T0/ctHb;
    CBV_w = CBV0*cbv_w;

    % Input
    tc  = x(1);
    tv  = x(2);
    Fpc = x(3);
    
    a = dynx(1);
    c = dynx(2);
    v = dynx(3);

    % Fixed parameters
    Sa = 0.98;
    pa = 0.3;
    alpha = 0.8;

    % Derived parameters
    pv = 1 - pa - Fpc;
    Sc  = Sa*(1/tc/alpha)*(1-exp(-alpha*tc)); % average capillary saturation
    Sv  = Sa*exp(-alpha*tc); % venous saturation

    % frequency response of RC Low pass filter (capillary)
    H_RC_LP_tmp     = (1./sqrt(1+(2*pi*f(1:n/2+1)*tc/exp(1)).^2)).*exp(-1i*atan(2*pi*f(1:n/2+1)*tc/exp(1)));
    H_RC_LP         = [H_RC_LP_tmp conj(H_RC_LP_tmp(end-1:-1:2))];
    % frequency response of Gaussian Low pass filter (venous)
    H_G_LP_tmp      = exp(-0.5*log(2)*(1.765*(tc+tv)*f(1:n/2+1)).^2).*exp(-1i*pi*(tc+tv)*f(1:n/2+1));
    H_G_LP          = [H_G_LP_tmp conj(H_G_LP_tmp(end-1:-1:2))];

    %% Dynamic cbf-cmro2

    CBVa_w = a*CBV_w;
    CBVc_w = c*CBV_w;
    CBVv_w = v*CBV_w;

    num = ((dO_w-dD_w)/T0)-(2*Sa-1)*CBVa_w/CBV0-(2*Sc-1)*CBVc_w/CBV0-(2*Sv-1)*CBVv_w/CBV0;
    denom = 2*((Sc/Sv)*(Sc-Sv)*Fpc*H_RC_LP+(Sa-Sv)*pv*H_G_LP);

    cbf_w = (conj(denom).*num.')./((conj(denom).*denom)+0);
    cbf_cmro2 = Fs*ifft(cbf_w,n,'symmetric');
    cbf_cmro2 = cbf_cmro2(1:L)'; 
    %% EOF
end

