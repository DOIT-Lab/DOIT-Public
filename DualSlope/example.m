%% Setup
% Giles Blaney (Giles.Blaney@tufts.edu) Fall 2019
clear; home;

%% Load and Parse Data
load('exampleData.mat');

% Example data from (Subject 4):
%   Blaney, G, Sassaroli, A, Pham, T, Fernandez, C, Fantini, S. Phase
%   dual?slopes in frequency?domain near?infrared spectroscopy for enhanced
%   sensitivity to brain tissue: First applications to human subjects. J.
%   Biophotonics. 2019;e201960018. https://doi.org/10.1002/jbio.201960018

% Physical locations: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                     ~   Src1     DetA   DetB     Src2   ~
%                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Data format: Det_Y(time, srcInd)
%              | srcInd | location | wavelength |
%              |      1 |     Src1 |          1 |
%              |      2 |     Src2 |          1 |
%              |      3 |     Src1 |          2 |
%              |      4 |     Src2 |          2 |

% Parsed format: Y(time, pair, wavelength)
%                Where "Y" is the datatype (I or phi)

I(:, :, 1)=[A_I(:, 1), B_I(:, 1), B_I(:, 2), A_I(:, 2)]; %counts
I(:, :, 2)=[A_I(:, 3), B_I(:, 3), B_I(:, 4), A_I(:, 4)]; %counts

phi(:, :, 1)=[A_phi(:, 1), B_phi(:, 1), B_phi(:, 2), A_phi(:, 2)]...
    *pi/180; %rad
phi(:, :, 2)=[A_phi(:, 3), B_phi(:, 3), B_phi(:, 4), A_phi(:, 4)]...
    *pi/180; %rad

%% Calculate O and D with Dual-Slope (DS)
opts.rho=[25, 35]; %mm
opts.nin=1.333;
opts.fmod=140.625e6; %Hz
opts.mua=0.01; %1/mm
opts.musp=1; %1/mm
opts.blInds=1:length(t);

initVar=NaN(length(t), size(I, 3));
dmua_DSI=initVar;
dmua_DSphi=initVar;
clear initVar;
for lInd=1:size(I, 3)
    dmua_DSI(:, lInd)=DSdmua(I(:, :, lInd), 'intensity', opts); %1/mm
    dmua_DSphi(:, lInd)=DSdmua(phi(:, :, lInd), 'phase', opts); %1/mm
end

[O_DSI, D_DSI]=mua2OandD(dmua_DSI*10, lambda); %uM
[O_DSphi, D_DSphi]=mua2OandD(dmua_DSphi*10, lambda); %uM

%% Plot
xl=[t(1), t(end)];
fs=1/median(diff(t));

figure(1); clf;
% Intensity
subplot(3, 2, 1);
lInd=1;
I0=mean(I(opts.blInds, :, lInd));
plot(t, lowpass((I(:, :, lInd)-I0)./I0+(-0.075:0.05:0.075), 0.2, fs));
xlim(xl);
legend('A1 (Short1)', 'B1 (Long1)', 'B2 (Short2)', 'A2 (Long2)',...
    'location', 'best');
ylabel('\DeltaI/I_0');
xlabel('t (sec)');
title(sprintf('Intensity Data %.0f nm', lambda(lInd)));

subplot(3, 2, 3);
lInd=2;
I0=mean(I(opts.blInds, :, lInd));
plot(t, lowpass((I(:, :, lInd)-I0)./I0+(-0.075:0.05:0.075), 0.2, fs));
xlim(xl);
ylabel('\DeltaI/I_0');
xlabel('t (sec)');
title(sprintf('Intensity Data %.0f nm', lambda(lInd)));

subplot(3, 2, 5);
plot(t, lowpass(O_DSI+0.3, 0.2, fs), '-r'); hold on;
plot(t, lowpass(D_DSI-0.3, 0.2, fs), '-b'); hold off;
xlim(xl);
legend('O', 'D', 'location', 'best');
ylabel('\Delta (\muM)');
xlabel('t (sec)');
title('Dual-Slope Intensity');

% Phase
subplot(3, 2, 2);
lInd=1;
phi0=phi(1, :, lInd);
plot(t, lowpass(phi(:, :, lInd)-phi0+(-0.015:0.01:0.015), 0.2, fs));
xlim(xl);
ylabel('\Delta\phi (rad)');
xlabel('t (sec)');
title(sprintf('Phase Data %.0f nm', lambda(lInd)));

subplot(3, 2, 4);
lInd=2;
phi0=phi(1, :, lInd);
plot(t, lowpass(phi(:, :, lInd)-phi0+(-0.015:0.01:0.015), 0.2, fs));
xlim(xl);
ylabel('\Delta\phi (rad)');
xlabel('t (sec)');
title(sprintf('Phase Data %.0f nm', lambda(lInd)));

subplot(3, 2, 6);
plot(t, lowpass(O_DSphi+0.3, 0.2, fs), '-r'); hold on;
plot(t, lowpass(D_DSphi-0.3, 0.2, fs), '-b'); hold off;
xlim(xl);
ylabel('\Delta (\muM)');
xlabel('t (sec)');
title('Dual-Slope Phase');