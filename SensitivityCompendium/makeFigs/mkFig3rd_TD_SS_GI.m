%% Setup
% Giles Blaney Ph.D. Fall 2023
clear; home;

svFig_params.Fontsize=10;
svFig_params.Linewidth=[];

%% Set Input Parameters
typNm='TD_SS_GI';
rs=[0, 0, 0]; %mm
rd=[25, 0, 0;
    35, 0, 0]; %mm

% Optical-props
optProp.nin= 1.333;
optProp.nout=1;
optProp.musp=1.1; %1/mm
optProp.g=   0.9;
optProp.mua= 0.011; %1/mm

% Name-values
NVA={
      'pert',     [10, 10, 2], ... %mm
        'dr',             0.5, ... %mm
        ...'dr',             0.1, ... %mm
        'xl',       [-10, 45], ... %mm
        'yl',       [-20, 20], ... %mm
        'zl',         [0, 30], ... %mm
      ...'fmod',           100e6, ... %Hz
       'ndt',            10e3, ...
      'tend',            10e3, ... %ps
        'tg',    [1500, 2000], ... %ps
       ...'tgE',    [   0, 1000], ... %ps
    'simTyp',            'DT', ...
     ...'detNA',             0.5, ...
        ...'np',             1e9, ...
    'usePar',            true, ...
    'FFTconv',           true, ...
    };

%% Calculate S
if ~exist([typNm '_3rd.mat'], 'file')
    addpath('../deps/');
    tic;
    [S, params, Svox]=makeS(typNm, rs, rd, optProp, NVA{:});
    runTime=toc;
    rmpath('../deps/');
    
    save([typNm '_3rd.mat']);
end

%% Plot 3rd Angle
load([typNm '_3rd.mat']);
Svox=reshape(filloutliers(Svox(:), 'linear',...
        'gesd', 'MaxNumOutliers', 1), size(Svox));
H=ones(NVA{2}/NVA{4});
plotData=convn(Svox, H, 'same');

plotPreset='none';

posSurf=0.04;
negSurf=-0.015;

svFig_params.doExport=true;
svFig_params.FigSz_cm=[17.75, NaN]; %cm

spcH=1.5; %cm
padH=0; 
spcV=1; %cm
padV=0.2; %cm
spcT=3.1; %cm

xl=NVA{6}; %mm
yl=NVA{8}; %mm
zl=[-5, NVA{10}(2)-5]; %mm
xsl=mean([rd(:, 1); rs(:, 1)]); %mm
ysl=0; %mm
zsl=10; %mm

addpath('plotHelpers/');
h=figure(310); clf;
h.Name=[typNm '_3rd'];
plot3rdAngle;

saveFig_senCompen;
makeLatexFigure;
rmpath('plotHelpers/');