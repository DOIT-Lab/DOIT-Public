%% Setup
% Giles Blaney Ph.D. Fall 2023
clear; home;

svFig_params.Fontsize=10;
svFig_params.Linewidth=[];

varNum=9;
drho_all=linspace(5, 45, varNum); %mm
mrho=27.5; %mm
varNm='drho';
titFormat='$\\Delta\\rho=%.0f$ mm';

%% Set Input Parameters
typNm='TD_SS_V';
rs=[0, 0, 0]; %mm
% rd=[35, 0, 0]; %mm

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
        'xl',       [-10, 60], ... %mm
        'yl',     [-5.5, 5.5], ... %mm
        'zl',         [0, 30], ... %mm
      'fmod',           100e6, ... %Hz
       'ndt',            10e3, ...
      'tend',            10e3, ... %ps
        ...'tg',    [1000, 2000], ... %ps
       ...'tgE',    [   0, 1000], ... %ps
    'simTyp',            'DT', ...
     ...'detNA',             0.5, ...
        ...'np',             1e9, ...
    'usePar',            true, ...
   'FFTconv',            true, ...
    };

%% Calculate S
Svox_all=NaN( ...
    round(diff(NVA{6})/NVA{4})+1, ...
    round(diff(NVA{8})/NVA{4})+1, ...
    round(diff(NVA{10})/NVA{4})+1, ...
    varNum);
rs_all=cell(varNum, 1);
rd_all=cell(varNum, 1);
for i=1:varNum
    rd=[...
        mrho-drho_all(i)/2, 0, 0;...
        mrho+drho_all(i)/2, 0, 0];
    
    if ~exist([typNm '_' varNm '_' num2str(i) '.mat'], 'file')
        addpath('../deps/');
        tic;
        [~, params, Svox]=makeS(typNm, rs, rd, optProp, NVA{:});
        runTime=toc;
        rmpath('../deps/');
        
        save([typNm '_' varNm '_' num2str(i) '.mat'], 'Svox', 'params',...
            'rs', 'rd');
    else
        load([typNm '_' varNm '_' num2str(i) '.mat']);
    end
    
    Svox_all(:, :, :, i)=Svox;
    rs_all{i}=rs;
    rd_all{i}=rd;
end

%% Plot
S_all=NaN( ...
    round(diff(NVA{6})/NVA{4})+1, ...
    round(diff(NVA{8})/NVA{4})+1, ...
    round(diff(NVA{10})/NVA{4})+1, ...
    varNum);
for i=1:varNum
    Svox=Svox_all(:, :, :, i);

    Svox=reshape(filloutliers(Svox(:), 'linear',...
            'gesd', 'MaxNumOutliers', 3), size(Svox));
    H=ones(NVA{2}/NVA{4});

    S_all(:, :, :, i)=convn(Svox, H, 'same');
end
plotData=S_all;

var_all=drho_all;

svFig_params.doExport=true;
svFig_params.FigSz_cm=[17.75, 20]; %cm

saNVA={'s', 0.03, 'mt', 0.14, 'mb', 0.055, 'ml', 0.07, 'mr', 0.03};

xl=NVA{6}; %mm
yl=NVA{8}; %mm
zl=[-5, NVA{10}(2)]; %mm
ysl=0; %mm

addpath('plotHelpers/');
h=figure(3122); clf;
h.Name=[typNm '_' varNm];
plot_varSP;

saveFig_senCompen;
makeLatexFigure_varSP;
rmpath('plotHelpers/');