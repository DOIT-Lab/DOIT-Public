%% Setup
clear; home;

rSrc=[5.5, 5, 0];
rDetA=[16.5, 5, 10];
rDetB=[28.5, 5, 10];
rDet=[rDetA; rDetB];

L=[45, 10, 10];

opts.omega=2*pi*100e6; %rad/sec
opts.BC='EBC';
opts.lMax=3;
opts.mMax=3;
opts.nMax=3;

mua=0.01;
musp=1;
mua0_ex=[0.005, 0.02];
musp0_ex=[0.5, 2];
ni=1.3;
no=1;

optProp.mua=mua;
optProp.musp=musp;
optProp.nin=ni;
optProp.nout=no;

IerrFrac=0.001;
Perr=0.1*pi/180;

kappa=1.2e-3;

%% Run Fits
T=NaN(size(rDet, 1), 1);
for i=1:size(rDet, 1)
    T(i)=Tcuv(rSrc, rDet(i, :), L, optProp, opts);
end
Terr=(abs(T)*IerrFrac)*exp(1i*Perr);

global xall;

optsCost=opts;
optsCost.nin=ni;
optsCost.nout=no;
optsCost.L=L;
optsCost.kappa=kappa;
muaRec=NaN(length(mua0_ex), length(musp0_ex));
muspRec=NaN(length(mua0_ex), length(musp0_ex));
costRec=NaN(length(mua0_ex), length(musp0_ex));
muaRecHist=cell(length(mua0_ex), length(musp0_ex));
muspRecHist=cell(length(mua0_ex), length(musp0_ex));
for mua0Ind=1:length(mua0_ex)
    for musp0Ind=1:length(musp0_ex)
        mu0=[mua0_ex(mua0Ind), musp0_ex(musp0Ind)];
        
        options=optimoptions('fmincon', 'OutputFcn', @outfun,...
            'Algorithm', 'interior-point');
%         options=optimoptions('lsqnonlin', 'OutputFcn', @outfun,...
%             'Algorithm', 'levenberg-marquardt');
        
        xall=[];
        [mu, cost]=fmincon(@cuvTcost,...
                mu0, [], [], [], [], [0, 0], [], [], options,...
                T, Terr, rSrc, rDet, optsCost);
%         [mu, cost]=lsqnonlin(@cuvTcost,...
%                 mu0, [0, 0], [], options,...
%                 T, Terr, rSrc, rDet, optsCost);
        
        muaRecHist{mua0Ind, musp0Ind}=xall(:, 1);
        muspRecHist{mua0Ind, musp0Ind}=xall(:, 2);
        muaRec(mua0Ind, musp0Ind)=mu(1);
        muspRec(mua0Ind, musp0Ind)=mu(2);
        costRec(mua0Ind, musp0Ind)=cost;
    end
end

%% Example Cost Space
mua_all=linspace(0.0025, 0.0225, 100);
musp_all=linspace(0.25, 2.25, 100);

costSP=NaN(length(mua_all), length(musp_all));
for i=1:length(mua_all)
    for j=1:length(musp_all)
        costSP(i, j)=...
            cuvTcost([mua_all(i), musp_all(j)],...
            T, Terr, rSrc, rDet, optsCost);
    end
end

%% Plot
% cl=[min(cost(:)), max(cost(:))];
cl=[1e-2, 1e4];
cTicks=[1e-1, 1e0, 1e1, 1e2, 1e3];

figure(90); clf;
subaxis(1, 1, 1, 'mr', 0.17, 'mb', 0.13, 'mt', 0.08);
pcolor(mua_all, musp_all, costSP.'); hold on;
contour(mua_all, musp_all, costSP.',...
    cTicks, 'w'); 
h1=plot(mua, musp, 'sk');
% h1=xline(mua, '--k');
% yline(musp, '--k');
h2=plot(NaN, NaN, ':ok');
cols=hsv(length(mua0_ex)*length(musp0_ex));
n=0;
for mua0Ind=1:length(mua0_ex)
    for musp0Ind=1:length(musp0_ex)
        n=n+1;
%         plot(muaRecHist{mua0Ind, musp0Ind}, muspRecHist{mua0Ind, musp0Ind},...
%             ':o', 'color', cols(n, :));
        plot(muaRecHist{mua0Ind, musp0Ind}, muspRecHist{mua0Ind, musp0Ind},...
            'o', 'color', cols(n, :));
        for i=1:(length(muaRecHist{mua0Ind, musp0Ind})-1)
            plot(muaRecHist{mua0Ind, musp0Ind}(i:(i+1)),...
                muspRecHist{mua0Ind, musp0Ind}(i:(i+1)),...
                ':', 'color', cols(n, :));
            
%             quiver(muaRecHist{mua0Ind, musp0Ind}(i),...
%                 muspRecHist{mua0Ind, musp0Ind}(i),...
%                 muaRecHist{mua0Ind, musp0Ind}(i+1)-muaRecHist{mua0Ind, musp0Ind}(i),...
%                 muspRecHist{mua0Ind, musp0Ind}(i+1)-muspRecHist{mua0Ind, musp0Ind}(i),...
%                 'color', cols(n, :), 'LineStyle', ':',...
%                 'AutoScale', 'on', 'AutoScaleFactor', 2)
        end
    end
end

hold off;
shading flat;
set(gca, 'ColorScale', 'log');
caxis(cl);
title(sprintf(...
    ['$\\mu_{a,true}=%.3f$ mm$^{-1}$, '...
    '$\\mu''_{s,true}=%.1f$ mm$^{-1}$'],...
    mua, musp),...
    'Interpreter', 'latex');

pos=get(gca, 'Position');
cb=colorbar;
set(gca, 'Position', pos);
set(cb, 'Ticks', cTicks);
ylabel(cb, '$\chi^2$', 'Interpreter', 'latex');

legend([h1, h2], 'True', 'Example fits iterations',...
    'Location', 'north');
ylabel('$\mu''_s$ (mm$^{-1}$)');
xlabel('$\mu_a$ (mm$^{-1}$)');


h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',6);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 8]]);
figName=['exampleFit'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');

%% Functions
function stop = outfun(x, optimvalues, state,...
    T, Terr, rSrc, rDet, optsCost)
    global xall;
    stop=false;
    if isequal(state, 'iter')
      xall=[xall; x];
    end
end