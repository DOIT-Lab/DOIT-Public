%% Setup
% Giles Blaney Ph.D. Fall 2023
clear; home;

svFig_params.Fontsize=10;
svFig_params.Linewidth=[];

%% Set Input Parameters
typNm='TD_SD_GI';
rs=[0, 0, 0]; %mm
rd=[0, 0, 0]; %mm

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
        'xl', [-15.25, 15.25], ... %mm
        'yl', [-15.25, 15.25], ... %mm
        'zl',         [0, 30], ... %mm
      ...'fmod',           100e6, ... %Hz
       'ndt',            1000, ...
      'tend',            2500, ... %ps
        'tg',    [1500, 2000], ... %ps
       ...'tgE',    [   0, 1000], ... %ps
    'simTyp',            'MC', ...
     'detNA',             0.5, ...
        'np',             1e9, ...
    ...'usePar',            true, ...
    ...'FFTconv',           true, ...
    };

%% Calculate S
if ~exist([typNm '_3rd_rho0_MC.mat'], 'file')
    addpath('../deps/');
    tic;
    [S, params, Svox]=makeS(typNm, rs, rd, optProp, NVA{:});
    runTime=toc;
    rmpath('../deps/');
    
    save([typNm '_3rd_rho0_MC.mat']);
end

%% Plot 3rd Angle Vox
load([typNm '_3rd_rho0_MC.mat']);
plotData=Svox;
plotPreset='vox';

posSurf=3e-5;
negSurf=NaN;

svFig_params.doExport=true;
svFig_params.FigSz_cm=[17.75, NaN]; %cm

spcH=1.5; %cm
padH=0; 
spcV=1.2; %cm
padV=0 %cm
spcT=3.5; %cm

xl=NVA{6}; %mm
yl=NVA{8}; %mm
zl=[-5, NVA{10}(2)-5]; %mm
xsl=mean([rd(:, 1); rs(:, 1)]); %mm
ysl=0; %mm
zsl=10; %mm

addpath('plotHelpers/');
h=figure(300); clf;
h.Name=[typNm '_3rd_rho0_MC'];
plot3rdAngle;

saveFig_senCompen;
makeLatexFigure;
rmpath('plotHelpers/');

%% Plot 3rd Angle
load([typNm '_3rd_rho0_MC.mat']);
plotData=S;
plotPreset='none';

posSurf=0.02;
negSurf=NaN;

svFig_params.doExport=true;
svFig_params.FigSz_cm=[17.75, NaN]; %cm

spcH=1.5; %cm
padH=0; 
spcV=1.5; %cm
padV=0; %cm
spcT=3.5; %cm

xl=NVA{6}; %mm
yl=NVA{8}; %mm
zl=[-5, NVA{10}(2)-5]; %mm
xsl=mean([rd(:, 1); rs(:, 1)]); %mm
ysl=0; %mm
zsl=10; %mm

addpath('plotHelpers/');
h=figure(301); clf;
h.Name=[typNm '_3rd_rho0_MC'];
plot3rdAngle;

saveFig_senCompen;
makeLatexFigure;
rmpath('plotHelpers/');