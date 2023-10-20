%% Setup
% Giles Blaney Ph.D. Fall 2023
clear; home;

svFig_params.Fontsize=10;
svFig_params.Linewidth=[];

%% Set Input Parameters
typNm='CW_SD_I';
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
        'dr',             0.1, ... %mm
        'xl', [-15.05, 15.05], ... %mm
        'yl', [-15.05, 15.05], ... %mm
        'zl',         [0, 30], ... %mm
      ...'fmod',           100e6, ... %Hz
       'ndt',              10, ...
      'tend',            10e3, ... %ps
        ...'tg',    [1000, 2000], ... %ps
       ...'tgE',    [   0, 1000], ... %ps
    'simTyp',            'MC', ...
     'detNA',             0.5, ...
        'np',             1e9, ...
    ...'usePar',            true, ...
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

%% Plot Voxelized 3rd Angle
load([typNm '_3rd_rho0_MC.mat']);
plotData=Svox;
plotPreset='vox';

posSurf=3e-5;
negSurf=NaN;

svFig_params.doExport=true;
svFig_params.FigSz_cm=[17.75, NaN]; %cm

spcH=1.5; %cm
padH=0; 
spcV=1; %cm
padV=0.5; %cm
spcT=3.5; %cm

xl=[-2, 2]; %mm
yl=[-2, 2]; %mm
zl=[-1, 2]; %mm
xsl=0; %mm
ysl=0; %mm
zsl=1; %mm

addpath('plotHelpers/');
h=figure(100); clf;
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
spcV=1; %cm
padV=0.5; %cm
spcT=3.2; %cm

xl=[-15, 15]; %mm
yl=[-15, 15]; %mm
zl=[-5, 10]; %mm
xsl=0; %mm
ysl=0; %mm
zsl=2; %mm

addpath('plotHelpers/');
h=figure(101); clf;
h.Name=[typNm '_3rd_rho0_MC'];
plot3rdAngle;

saveFig_senCompen;
makeLatexFigure;
rmpath('plotHelpers/');