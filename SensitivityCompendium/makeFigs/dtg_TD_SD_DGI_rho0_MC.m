%% Setup
% Giles Blaney Ph.D. Fall 2023
clear; home;

svFig_params.Fontsize=10;
svFig_params.Linewidth=[];

mtg=1250; %ps
dt=500; %ps
dtg_all=500:300:1100; %ps
varNum=length(dtg_all);
varNm='dtg';
titFormat='$\\Delta (\\bar{t})_{gates}=%.0f$ ps';

%% Set Input Parameters
typNm='TD_SD_DGI';
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
         'xl', [-15.25, 15.25], ... %mm
         'yl', [-15.25, 15.25], ... %mm
         'zl',         [0, 30], ... %mm
       ...'fmod',           100e6, ... %Hz
        'ndt',            2080*4, ...
       'tend',            5200, ... %ps
           'tg',    [NaN, NaN], ... %ps
          'tgE',    [NaN, NaN], ... %ps
     'simTyp',            'MC', ...
      'detNA',             0.5, ...
         'np',             1e9, ...
      ...'usePar',           true, ...
     ...'FFTconv',           true,...
    };

%% Calculate S
S_all=NaN( ...
    round(diff(NVA{6})/NVA{4}), ...
    round(diff(NVA{8})/NVA{4}), ...
    round(diff(NVA{10})/NVA{4}), ...
    varNum);
rs_all=cell(varNum, 1);
rd_all=cell(varNum, 1);
for i=1:varNum
    mtEarly=mtg-dtg_all(i)/2;
    mtLate=mtg+dtg_all(i)/2;
    
    NVA{18}=[mtEarly-dt/2, mtEarly+dt/2];
    NVA{16}=[mtLate-dt/2, mtLate+dt/2];

    if ~exist([typNm '_' varNm '_' num2str(i) '_rho0_MC.mat'], 'file')
        addpath('../deps/');
        tic;
        [S, params, ~]=makeS(typNm, rs, rd, optProp, NVA{:});
        runTime=toc;
        rmpath('../deps/');
        
        save([typNm '_' varNm '_' num2str(i) '_rho0_MC.mat'], 'S', 'params',...
            'rs', 'rd');
    else
        load([typNm '_' varNm '_' num2str(i) '_rho0_MC.mat']);
    end
    
    S_all(:, :, :, i)=S;
    rs_all{i}=rs;
    rd_all{i}=rd;
end

%% Plot
plotData=S_all;
var_all=dtg_all;

svFig_params.doExport=true;
svFig_params.FigSz_cm=[17.75, 10]; %cm

saNVA={'s', 0.03, 'mt', 0.17, 'mb', 0.055, 'ml', 0.07, 'mr', 0.03};

xl=NVA{6}; %mm
yl=NVA{8}; %mm
zl=[-5, NVA{10}(2)]; %mm
ysl=0; %mm

addpath('plotHelpers/');
h=figure(3037); clf;
h.Name=[typNm '_' varNm '_rho0_MC'];
plot_varSP;

saveFig_senCompen;
makeLatexFigure_varSP;
rmpath('plotHelpers/');