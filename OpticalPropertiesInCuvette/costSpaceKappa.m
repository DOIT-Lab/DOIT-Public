%% Setup
clear; home;

rSrc=[5.5, 5, 0];
rDetA=[16.5, 5, 10];
rDetB=[28.5, 5, 10];

L=[45, 10, 10];

opts.omega=2*pi*100e6; %rad/sec
opts.BC='EBC';
opts.lMax=3;
opts.mMax=3;
opts.nMax=3;

mua=0.01;
musp=1;
ni=1.3;
no=1;

optProp.mua=mua;
optProp.musp=musp;
optProp.nin=ni;
optProp.nout=no;

IerrFrac=0.001;
Perr=0.1*pi/180;

%% Sim T
rDet=[rDetA; rDetB];
T=NaN(size(rDet, 1), 1);
for i=1:size(rDet, 1)
    T(i)=Tcuv(rSrc, rDet(i, :), L, optProp, opts);
end
Terr=(abs(T)*IerrFrac)*exp(1i*Perr);

%% Find Kappa
mua_all=linspace(0.5*mua, 1.5*mua, 50);
musp_all=linspace(0.5*musp, 1.5*musp, 50);

optsKapCost=opts;
optsKapCost.mua_all=mua_all;
optsKapCost.musp_all=musp_all;
optsKapCost.nin=ni;
optsKapCost.nout=no;
optsKapCost.L=L;
optsKapCost.Tmeas=T;
optsKapCost.Terr=Terr;
optsKapCost.rSrc=rSrc;
optsKapCost.rDet=rDet;

kappaLo=0.00003;
kappaUp=0.03;
kappa0=0.001;
[kappa, kappaC]=fmincon(@kappaCost,...
    kappa0, [], [], [], [], kappaLo, kappaUp, [], [], optsKapCost);

%% Example Cost Space
mua_all=linspace(0.5*mua, 1.5*mua, 100);
musp_all=linspace(0.5*musp, 1.5*musp, 100);

optsCost=opts;
optsCost.nin=ni;
optsCost.nout=no;
optsCost.L=L;

optsCost.kappa=kappaLo;
costLo=NaN(length(mua_all), length(musp_all));
for i=1:length(mua_all)
    for j=1:length(musp_all)
        costLo(i, j)=cuvTcost([mua_all(i), musp_all(j)],...
            T, Terr, rSrc, rDet, optsCost);
    end
end
kappaClo=kappaCost(kappaLo, optsKapCost);

optsCost.kappa=kappaUp;
costUp=NaN(length(mua_all), length(musp_all));
for i=1:length(mua_all)
    for j=1:length(musp_all)
        costUp(i, j)=cuvTcost([mua_all(i), musp_all(j)],...
            T, Terr, rSrc, rDet, optsCost);
    end
end
kappaCup=kappaCost(kappaUp, optsKapCost);

optsCost.kappa=kappa;
cost=NaN(length(mua_all), length(musp_all));
for i=1:length(mua_all)
    for j=1:length(musp_all)
        cost(i, j)=cuvTcost([mua_all(i), musp_all(j)],...
            T, Terr, rSrc, rDet, optsCost);
    end
end

%% Plot
cl=[0, 1];

figure(70); clf;
subaxis(1, 3, 1, 'mr', 0.15, 'mt', 0.15);
pcolor(mua_all/mua, musp_all/musp, costLo.'/max(costLo(:))); hold on;
contour(mua_all/mua, musp_all/musp,...
    costLo.'/max(costLo(:)),...
    [1, 1]*quantile(costLo(:).'/max(costLo(:)), 0.05), 'w');
shading flat; hold off; axis equal tight;
caxis(cl);
ylabel('$\mu''_s/\mu''_{s,true}$');
xlabel('$\mu_a/\mu_{a,true}$');
tmpStr=[strrep(sprintf('%.1e', kappaLo),...
        'e', '\times 10^{'), '}'];
title(sprintf('\\textbf{(a)} $\\kappa=%s$\n$P^2/A=%.0f$',...
    tmpStr, kappaClo),...
    'Interpreter', 'latex');

subaxis(1, 3, 2);
pcolor(mua_all/mua, musp_all/musp, cost.'/max(cost(:))); hold on;
contour(mua_all/mua, musp_all/musp,...
    cost.'/max(cost(:)),...
    [1, 1]*quantile(cost(:).'/max(cost(:)), 0.05), 'w');
h1=plot(NaN, NaN, '-k');
shading flat; hold off; axis equal tight;
caxis(cl);
set(gca, 'YTickLabel', {});
xlabel('$\mu_a/\mu_{a,true}$');
tmpStr=[strrep(sprintf('%.1e', kappa),...
        'e', '\times 10^{'), '}'];
leg=legend(h1, '$\chi^2_{image}$ 5$^{th}$ quantile',...
    'Location', 'northwest', 'Interpreter', 'latex');
set(leg, 'Position', leg.Position+[-0.15, 0, 0, 0]);
title(sprintf('\\textbf{(b)} $\\kappa=%s$\n$P^2/A=%.0f$',...
    tmpStr, kappaC),...
    'Interpreter', 'latex');

subaxis(1, 3, 3);
pcolor(mua_all/mua, musp_all/musp, costUp.'/max(costUp(:))); hold on;
contour(mua_all/mua, musp_all/musp,...
    costUp.'/max(costUp(:)),...
    [1, 1]*quantile(costUp(:).'/max(costUp(:)), 0.05), 'w');
pos=get(gca, 'Position');
cb=colorbar;
set(gca, 'Position', pos);
shading flat; hold off; axis equal tight;
caxis(cl);
set(gca, 'YTickLabel', {});
xlabel('$\mu_a/\mu_{a,true}$');
ylabel(cb, '$\chi^2/\max(\chi^2_{image})$',...
    'Interpreter', 'latex');
tmpStr=[strrep(sprintf('%.1e', kappaUp),...
        'e', '\times 10^{'), '}'];
title(sprintf('\\textbf{(c)} $\\kappa=%s$\n$P^2/A=%.0f$',...
    tmpStr, kappaCup),...
    'Interpreter', 'latex');

sgtitle(sprintf(...
    ['$\\mu_{a,true}=%.3f$ mm$^{-1}$, '...
    '$\\mu''_{s,true}=%.1f$ mm$^{-1}$'],...
    mua, musp),...
    'Interpreter', 'latex');


h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',6);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 6]]);
figName=['exKappa'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');

%% Functions
function kappaCostVal = kappaCost(kappa, opts)

    opts.kappa=kappa;
    
    mua_all=opts.mua_all;
    musp_all=opts.musp_all;
    Tmeas=opts.Tmeas;
    Terr=opts.Terr;
    rSrc=opts.rSrc;
    rDet=opts.rDet;
    
    cost=NaN(length(opts.mua_all), length(opts.musp_all));
    n_mua=length(opts.mua_all);
    n_musp=length(opts.musp_all);
    parfor i=1:n_mua
        for j=1:n_musp
            cost(i, j)=cuvTcost([mua_all(i), musp_all(j)],...
                Tmeas, Terr, rSrc, rDet, opts);
        end
    end
    costNorm=cost/max(cost(:));
    
    [II, JJ]=meshgrid(linspace(0, 1, length(mua_all)),...
        linspace(0, 1, length(musp_all)));
    [III, JJJ]=meshgrid(...
        linspace(0, 1, 1000),...
        linspace(0, 1, 1000));
    cost_up=interp2(II, JJ, costNorm, III, JJJ, 'spline');

    figure(999); clf;
    pcolor(III(1, :), JJJ(:, 1), cost_up.'); hold on;
    M=contour(III(1, :), JJJ(:, 1), cost_up.',...
        [1, 1]*quantile(cost_up(:), 0.05), 'w'); hold off; shading flat;
    axis equal tight off;
    shp=alphaShape(M(1, 2:M(2, 1)).', M(2, 2:M(2, 1)).');

    kappaCostVal=perimeter(shp)^2/area(shp);

    title(sprintf('$\\kappa=%f$\n$P^2/A=%f$', kappa, kappaCostVal),...
        'Interpreter', 'latex');
    drawnow;
end