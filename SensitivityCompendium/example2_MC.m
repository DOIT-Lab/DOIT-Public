%% Setup
% Giles Blaney Ph.D. Summer 2023
clear; home;

addpath('deps/');

% % Single-Distance Examples
typNm='CW_SD_I';
% typNm='TD_SD_GI';
% typNm='TD_SD_DGI';
rho1=0; %mm
rs=[0, 0, 0];
rd=[rho1, 0, 0];

% % % Single-Slope Examples
% typNm='CW_SS_I';
% % typNm='TD_SS_GI';
% rho1=0; %mm
% rho2=5; %mm
% rs=[0, 0, 0];
% rd=[rho1, 0, 0; rho2, 0, 0];

% % % Dual-Slope Examples
% typNm='CW_DS_I';
% % typNm='TD_DS_GI';
% rho1=0; %mm
% rho2=5; %mm
% rs=[0, 0, 0; rho1+rho2, 0, 0];
% rd=[rho1, 0, 0; rho2, 0, 0];

%% Set Input Parameters
pert=[1, 1, 1]; %mm
dr=0.5; %mm
xl=[-15.5, 20.5];
yl=[-15.5, 15.5];
zl=[0, 15];

optProp.nin=1.333;
optProp.nout=1;
optProp.musp=1.1; %1/mm
optProp.g=0.9;
optProp.mua=0.011; %1/mm

% For Gated TD
tg=[1000, 2000]; %ps

%% Calculate S
[S, params]=makeS(typNm, rs, rd, optProp,...
    'pert', pert, 'dr', dr, 'xl', xl, 'yl', yl', 'zl', zl, ...
    'tg', tg, 'simTyp', 'MC');

%% Plot
[Splane, pp]=sliceS(S, params, 'y', 0);

[cl, cols]=makeCL(S, 'qants', [0.01, 0.99]);

h=figure(20); clf;
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