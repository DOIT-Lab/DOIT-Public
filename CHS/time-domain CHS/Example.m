close all; clear; clc;

%% Setup

% Load sasmple of absolute totalhemoglobin time traces
load('testdata');
[~,t0idx] = min(abs(t-0)); % t = 0s is at the step release
T0 = T(t0idx); % microM
dT = T-T0; % microM
L = length(dT);

% Model parameters
ctHb = 2300; % microM
tc = 0.7; % s
tv = 5; % s
Fpc = 0.4;
fc = 0.05; % Hz
k = 3;
x = [tc tv Fpc fc k];
dynx = [0.2 0 0.8];

%% Forward model
[O,D] = forward_TD(x,dynx,T,T0,ctHb,Fs);

figure(1), 
subplot(1,3,1),hold on
plot(t,O,'r');
axis tight;
ylim(round([0.99*min(O) 1.01*max(O)],2));
plot([0 0],ylim,'--k');
xlabel('Time (s)');
ylabel('O (\muM)');
title('Oxy-hemoglobin conc.');
subplot(1,3,2), hold on
plot(t,D,'b');
axis tight;
ylim(round([0.9*min(D) 1.1*max(D)],2));
plot([0 0],ylim,'--k');
xlabel('Time (s)');
ylabel('D (\muM)');
title('Deoxy-hemoglobin conc.');
subplot(1,3,3), hold on
plot(t,T,'g');
axis tight;
ylim(round([0.99*min(T) 1.01*max(T)],2));
plot([0 0],ylim,'--k');
xlabel('Time (s)');
ylabel('T (\muM)');
title('Total-hemoglobin conc.');
sgtitle('Forward time-domain CHS model');

%% Inverse model

% Add noise to the data
O = O+0.4*randn(L,1); % microM
D = D+0.04*randn(L,1); % microM
O0 = mean(O(t<=0)); % microM
D0 = mean(D(t<=0)); % microM
T0 = O0+D0;

% Low-pass filtered
b = LPfilt2(Fs,0.2,0.5);
dOf = filtfilt(b,1,O-O0);
dDf = filtfilt(b,1,D-D0);
Of = dOf+O0;
Df = dDf+D0;

% Using temporal O & D during 10 s after the step release for the fit
timestep = t((t>=0)&(t<=10));
Ostep = Of((t>=0)&(t<=10));
Dstep = Df((t>=0)&(t<=10));
Ostep_err = std(O((t>=0)&(t<=10))-Ostep);
Dstep_err = std(D((t>=0)&(t<=10))-Dstep);
data = [Ostep,Dstep];
data_err = [Ostep_err,Dstep_err];
[xfit,xfit_s,datafit,chi2r] = Inverse_TD(dynx,ctHb,data,[1,1],Fs);

figure(2),
subplot(1,2,1), hold on
plot(timestep,Ostep,'*r');
plot(timestep,datafit(:,1),'k');
legend('data','fit');
xlabel('Time (s)');
ylabel('O (\muM)');
title('Oxy-hemoglobin conc.');
subplot(1,2,2), hold on
plot(timestep,Dstep,'*b');
plot(timestep,datafit(:,2),'k');
legend('data','fit');
xlabel('Time (s)');
ylabel('D (\muM)');
title('Deoxy-hemoglobin conc.');
sgtitle('Inverse time-domain CHS model');

parlab = {'t^{(c)} (s)','t^{(v)} (s)','F^{(c)}CBV_0^{(c)}/CBV_0','f^{(c)} (Hz)','k'};
parlim = [[0 2];[0 7];[0 1];[0 0.2];[1 5]];
figure(3),
for i = 1:length(xfit)
    subplot(1,5,i), hold on
    plot([0 2],[x(i) x(i)],'k');
    errorbar(1,xfit(i),xfit_s(i),'o');
    axis tight;
    ylim(parlim(i,:));
    set(gca,'XTick',{});
    if i==length(xfit)
        legend('True','Fitted');
    end
    title(parlab{i},'fontweight','normal');
end
sgtitle('Fitted parameters');


