%% Setup
% Giles Blaney Ph.D. Summer 2023
clear; home;

addpath('deps/');

% % Single-Distance Examples
% typNm='CW_SD_I';
% typNm='FD_SD_I';
% typNm='FD_SD_P';
typNm='TD_SD_GI';
% typNm='TD_SD_DGI';
% typNm='TD_SD_T';
% typNm='TD_SD_V';
rho2=35; %mm
rs=[0, 0, 0];
rd=[rho2, 0, 0];

% % % Single-Slope Examples
% typNm='CW_SS_I';
% % typNm='FD_SS_I';
% % typNm='FD_SS_P';
% % typNm='TD_SS_GI';
% % typNm='TD_SS_T';
% % typNm='TD_SS_V';
% rho1=25; %mm
% rho2=35; %mm
% rs=[0, 0, 0];
% rd=[rho1, 0, 0; rho2, 0, 0];

% % % Dual-Slope Examples
% typNm='CW_DS_I';
% % typNm='FD_DS_I';
% % typNm='FD_DS_P';
% % typNm='TD_DS_GI';
% % typNm='TD_DS_T';
% % typNm='TD_DS_V';
% rho1=25; %mm
% rho2=35; %mm
% rs=[0, 0, 0; rho1+rho2, 0, 0];
% rd=[rho1, 0, 0; rho2, 0, 0];

%% Set Input Parameters
pert=[1, 1, 1]; %mm
dr=1; %mm
xl=[-10, 70];
yl=[0, 0];
zl=[0, 25];

optProp.nin=1.333;
optProp.nout=1;
optProp.musp=1.1; %1/mm
optProp.mua=0.011; %1/mm

% For FD
fmod=100e6; %Hz

% For Gated TD
tg=[1000, 2000]; %ps

%% Calculate S
[S, params]=makeS(typNm, rs, rd, optProp,...
    'pert', pert, 'dr', dr, 'xl', xl, 'yl', yl', 'zl', zl, ...
    'fmod', fmod, 'tg', tg, ...
    'usePar', false, 'FFTconv', true);

%% Plot
[Splane, pp]=sliceS(S, params, 'y', 0);

[cl, cols]=makeCL(S);

h=figure(99); clf;
h.Name=typNm;
imagesc(pp.horzAx, pp.vertAx, Splane); hold on;
clim(cl);
set(gca, 'Colormap', cols);
axis equal tight;
cb=colorbar;
ylabel(cb, '$S$', 'Interpreter', 'latex');
contour(pp.horzAx, pp.vertAx, Splane, cb.Ticks, '--',...
    'Color', [1, 1, 1]*0.5);  hold off;
xlabel(pp.horzNm, 'Interpreter', 'latex');
ylabel(pp.vertNm, 'Interpreter', 'latex');
title(typNm, 'Interpreter', 'none');

%% Cleanup
rmpath('deps/');