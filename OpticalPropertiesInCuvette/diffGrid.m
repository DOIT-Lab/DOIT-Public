%% Setup 
clear; home;

rSrc=[5.5, 5, 0];
rDetA=[16.5, 5, 10];
rDetB=[28.5, 5, 10];

L=[45, 10, 10];

mua=0.01;
musp=1;
mua_all=linspace(0, 0.05, 100);
musp_all=linspace(0.5, 5, 100);

ni=1.3;
no=1;
ni_all=linspace(1, 2, 100);
no_all=linspace(1, 2, 100);

opts.omega=2*pi*100e6; %rad/sec
opts.BC='EBC';
opts.lMax=3;
opts.mMax=3;
opts.nMax=3;

%% Run T
TA=NaN(length(mua_all), length(musp_all));
TB=NaN(length(mua_all), length(musp_all));
for muaInd=1:length(mua_all)
    for muspInd=1:length(musp_all)
        optProp.mua=mua_all(muaInd);
        optProp.musp=musp_all(muspInd);
        optProp.nin=ni;
        optProp.nout=no;
        TA(muaInd, muspInd)=Tcuv(rSrc, rDetA, L,...
            optProp, opts);
        TB(muaInd, muspInd)=Tcuv(rSrc, rDetB, L,...
            optProp, opts);
    end
end

Idiff=log(abs(TB))-log(abs(TA));
Pdiff=angle(TB)-angle(TA);


TAn=NaN(length(ni_all), length(no_all));
TBn=NaN(length(ni_all), length(no_all));
for niInd=1:length(ni_all)
    for noInd=1:length(no_all)
        optProp.mua=mua;
        optProp.musp=musp;
        optProp.nin=ni_all(niInd);
        optProp.nout=no_all(noInd);
        TAn(niInd, noInd)=Tcuv(rSrc, rDetA, L,...
            optProp, opts);
        TBn(niInd, noInd)=Tcuv(rSrc, rDetB, L,...
            optProp, opts);
    end
end

Indiff=log(abs(TBn))-log(abs(TAn));
Pndiff=angle(TBn)-angle(TAn);

%% Plot
Iticks=-10:2:-4;
Pticks=0.1:0.1:0.5;
clI=[min([Idiff(:); Indiff(:)]), max([Idiff(:); Indiff(:)])];
clP=[min([Pdiff(:); Pndiff(:)]), max([Pdiff(:); Pndiff(:)])];

figure(40); clf;

subaxis(2, 2, 2, 'sv', 0.07, 'sh', 0.1,...
    'ml', 0.08, 'mr', 0.15, 'mt', 0.1);
pcolor(ni_all, no_all, Indiff.'); shading flat; hold on;
caxis(clI);
pos=get(gca, 'Position');
cb1=colorbar;
set(gca, 'Position', pos);
cb1.Ticks=Iticks;
contour(ni_all, no_all, Indiff.', Iticks, 'w');
text(1.05, 1.9, '\textbf{(b)}', 'Interpreter', 'latex'); hold off;
set(gca, 'YDir', 'normal');
set(gca, 'XTickLabel', {});
ylabel('$n_o$');
title(sprintf('$\\mu_a=%.3f$ mm$^{-1}$, $\\mu''_s=%.2f$ mm$^{-1}$\nfor \\textbf{(b)} and \\textbf{(d)}',...
    mua, musp), 'Interpreter', 'latex');
ylabel(cb1, '$\ln|$SR$\{\widetilde{T}\}|_{\mathrm{\texttt{1AB}}}=\ln|\widetilde{T}_{\mathrm{\texttt{1B}}}|-\ln|\widetilde{T}_{\mathrm{\texttt{1A}}}|$',...
    'Interpreter', 'latex');
% ylabel(cb1, '$\ln|\widetilde{T}_{1B}|-\ln|\widetilde{T}_{1A}|$',...
%     'Interpreter', 'latex');

subaxis(2, 2, 4);
pcolor(ni_all, no_all, Pndiff.'); shading flat; hold on;
caxis(clP);
pos=get(gca, 'Position');
cb2=colorbar;
set(gca, 'Position', pos);
cb2.Ticks=Pticks;
contour(ni_all, no_all, Pndiff.', Pticks, 'w');
text(1.05, 1.9, '\textbf{(d)}', 'Interpreter', 'latex'); hold off;
set(gca, 'YDir', 'normal');
ylabel('$n_o$');
xlabel('$n_i$');
ylabel(cb2, '$\angle$SR$\{\widetilde{T}\}_{\mathrm{\texttt{1AB}}}=\angle\widetilde{T}_{\mathrm{\texttt{1B}}}-\angle\widetilde{T}_{\mathrm{\texttt{1A}}}$ (rad)',...
    'Interpreter', 'latex');


subaxis(2, 2, 1);
pcolor(mua_all, musp_all, Idiff.'); shading flat; hold on;
caxis(clI);
contour(mua_all, musp_all, Idiff.', Iticks, 'w'); 
text(0.0015, 4.6, '\textbf{(a)}', 'Interpreter', 'latex'); hold off;
set(gca, 'YDir', 'normal');
set(gca, 'XTickLabel', {});
ylabel('$\mu''_s$ (mm$^{-1}$)');
title(sprintf('$n_i=%.1f$, $n_o=%.1f$\nfor \\textbf{(a)} and \\textbf{(c)}',...
    ni, no), 'Interpreter', 'latex');

subaxis(2, 2, 3);
pcolor(mua_all, musp_all, Pdiff.'); shading flat; hold on;
caxis(clP);
contour(mua_all, musp_all, Pdiff.', Pticks, 'w'); 
text(0.0015, 4.6, '\textbf{(c)}', 'Interpreter', 'latex'); hold off;
set(gca, 'YDir', 'normal');
xlabel('$\mu_a$ (mm$^{-1}$)');
ylabel('$\mu''_s$ (mm$^{-1}$)');

h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',5);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 13.5]]);
figName=['diffGrd'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');