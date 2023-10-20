%% Setup
% Giles Blaney Ph.D. Fall 2023
clear; home;

svFig_params.Fontsize=10;
svFig_params.Linewidth=[];

varNum=8;
mrho_all=linspace(10, 45, varNum); %mm
drho=10; %mm
varNm='mrho';
titFormat='$\\bar{\\rho}=%.0f$ mm';

%% Set Input Parameters
typNm='CW_DS_I';
% rs=[-mrho, 0, 0;...
%     mrho, 0, 0]; %mm
rd=[...
    -drho/2, 0, 0;...
    drho/2, 0, 0];

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
        'xl',       [-max(mrho_all)-10, max(mrho_all)+10], ... %mm
        'yl',       [-5.1, 5.1], ... %mm
        'zl',         [0, 30], ... %mm
      ...'fmod',           100e6, ... %Hz
       ...'ndt',              10, ...
      ...'tend',            10e3, ... %ps
        ...'tg',    [1000, 2000], ... %ps
       ...'tgE',    [   0, 1000], ... %ps
    'simTyp',            'DT', ...
     ...'detNA',             0.5, ...
        ...'np',             1e9, ...
    ...'usePar',            true, ...
    };

%% Calculate S
S_all=NaN( ...
    round(diff(NVA{6})/NVA{4})+1, ...
    round(diff(NVA{8})/NVA{4})+1, ...
    round(diff(NVA{10})/NVA{4})+1, ...
    varNum);
rs_all=cell(varNum, 1);
rd_all=cell(varNum, 1);
for i=1:varNum
    rs=[...
        -mrho_all(i), 0, 0;...
        mrho_all(i), 0, 0];
    
    if ~exist([typNm '_' varNm '_' num2str(i) '.mat'], 'file')
        addpath('../deps/');
        tic;
        [S, params, ~]=makeS(typNm, rs, rd, optProp, NVA{:});
        runTime=toc;
        rmpath('../deps/');
        
        save([typNm '_' varNm '_' num2str(i) '.mat'], 'S', 'params',...
            'rs', 'rd');
    else
        load([typNm '_' varNm '_' num2str(i) '.mat']);
    end
    
    S_all(:, :, :, i)=S;
    rs_all{i}=rs;
    rd_all{i}=rd;
end

%% Plot
plotData=S_all;
var_all=mrho_all;

svFig_params.doExport=true;
svFig_params.FigSz_cm=[17.75, 17]; %cm

saNVA={'s', 0.03, 'mt', 0.14, 'mb', 0.055, 'ml', 0.07, 'mr', 0.03};

xl=NVA{6}; %mm
yl=NVA{8}; %mm
zl=[-5, NVA{10}(2)]; %mm
ysl=0; %mm

addpath('plotHelpers/');
h=figure(1151); clf;
h.Name=[typNm '_' varNm];
plot_varSP;

saveFig_senCompen;
makeLatexFigure_varSP;
rmpath('plotHelpers/');