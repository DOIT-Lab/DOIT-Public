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

%% Run T
y=0:0.1:L(2);
x=0:0.1:L(1);

optProp.mua=0.01; %1/mm
optProp.musp=1; %1/mm

T=NaN(length(y), length(x));
for yInd=1:length(y)
    for xInd=1:length(x)
        rDet=[x(xInd), y(yInd)];
        T(yInd, xInd)=Tcuv(rSrc, rDet, L, optProp, opts);
    end
end

%% Plot
cols=lines(7);

ampLim=[1e-9, 1e-3];
phLim=[0.1, 0.7];

[~, y0Ind]=min(abs(y-mean(y)));

figure(20); clf;
subaxis(4, 1, 1, 's', 0.01,...
    'mt', 0.06, 'ml', 0.11, 'mr', 0.1);
imagesc(x-mean(x), y-mean(y), abs(T)); hold on;
plot(rSrc(1)-mean(x), rSrc(2)-mean(y), '+w'); 
plot(rSrc(1)-mean(x), rSrc(2)-mean(y), '.k'); 
plot(rDetA(1)-mean(x), rDetA(2)-mean(y), 'o', 'color', cols(1, :));
plot(rDetA(1)-mean(x), rDetA(2)-mean(y), '.k');
plot(rDetB(1)-mean(x), rDetB(2)-mean(y), 'o', 'color', cols(2, :));
plot(rDetB(1)-mean(x), rDetB(2)-mean(y), '.k');
text(-22, 3, '\textbf{(a)}', 'Interpreter', 'latex');
hold off;
axis equal tight;
caxis(ampLim);
cb=colorbar;
cb.Ticks=[1e-9, 1e-6, 1e-3];
ax1=gca;
ax1.Colormap=turbo;
ax1.XTickLabel={};
ax1.ColorScale='log';
ax1.YDir='normal';
for i=1:length(ax1.XTickLabel)
    ax1.XTickLabel{i}=['$' ax1.XTickLabel{i} '$'];
    ax1.XTick(i)=ax1.XTick(i);
end
for i=1:length(ax1.YTickLabel)
    ax1.YTickLabel{i}=['$' ax1.YTickLabel{i} '$'];
    ax1.YTick(i)=ax1.YTick(i);
end
ylabel('y (mm)');
ylabel(cb, '$|\widetilde{T}|$ (mm$^{-2}$)', 'Interpreter', 'latex');
title(sprintf('$z=10$ mm'),...
    'Interpreter', 'latex');

subaxis(4, 1, 2, 's', 0.01);
imagesc(x-mean(x), y-mean(y), angle(T)); hold on;
plot(rSrc(1)-mean(x), rSrc(2)-mean(y), '+w'); 
plot(rSrc(1)-mean(x), rSrc(2)-mean(y), '.k'); 
plot(rDetA(1)-mean(x), rDetA(2)-mean(y), 'o', 'color', cols(1, :));
plot(rDetA(1)-mean(x), rDetA(2)-mean(y), '.k');
plot(rDetB(1)-mean(x), rDetB(2)-mean(y), 'o', 'color', cols(2, :));
plot(rDetB(1)-mean(x), rDetB(2)-mean(y), '.k');
text(-22, 3, '\textbf{(b)}', 'Interpreter', 'latex');
hold off;
axis equal tight;
caxis(phLim);
cb=colorbar;
cb.Ticks=[0.2, 0.4, 0.6];
ax2=gca;
ax2.Colormap=parula;
ax2.YDir='normal';
for i=1:length(ax2.XTickLabel)
    ax2.XTickLabel{i}=['$' ax2.XTickLabel{i} '$'];
    ax2.XTick(i)=ax2.XTick(i);
end
for i=1:length(ax2.YTickLabel)
    ax2.YTickLabel{i}=['$' ax2.YTickLabel{i} '$'];
    ax2.YTick(i)=ax2.YTick(i);
end
xlabel('x (mm)');
ylabel('y (mm)');
ylabel(cb, '$\angle\widetilde{T}$ (rad)', 'Interpreter', 'latex');

% pTmp=ax1.Position;
% ax2.Position(3)=pTmp(3);
% ax1.Position(3)=pTmp(3);
% ax2.Position(3)=pTmp(3);

subaxis(4, 1, [3, 4], 's', 0.25);
yyaxis left;
semilogy(x-mean(x), abs(T(y0Ind, :)), '-b'); hold on;
h2=xline(rSrc(1)-mean(x), ':k');
h4=xline(rDetA(1)-mean(x), ':', 'color', cols(1, :));
h6=xline(rDetB(1)-mean(x), ':', 'color', cols(2, :));
h1=plot(NaN, NaN, '+k');
h3=plot(NaN, NaN, 'o', 'color', cols(1, :));
h5=plot(NaN, NaN, 'o', 'color', cols(2, :));
text(-22, 5e-5, '\textbf{(c)}', 'Interpreter', 'latex');
hold off;
ylim(ampLim);
xlim([x(1), x(end)]-mean(x));
ax=gca;
ax.YTick=[1e-9, 1e-6, 1e-3];
for i=1:length(ax.XTickLabel)
    ax.XTickLabel{i}=['$' ax.XTickLabel{i} '$'];
    ax.XTick(i)=ax.XTick(i);
end
for i=1:length(ax.YTickLabel)
    ax.YTickLabel{i}=['$' ax.YTickLabel{i} '$'];
    ax.YTick(i)=ax.YTick(i);
end
ylabel('$|\widetilde{T}|$ (mm$^{-2}$)', 'Interpreter', 'latex');

yyaxis right;
plot(x-mean(x), angle(T(y0Ind, :)), '-r');
ylim([phLim(1), phLim(2)*1.5]);
axr=gca;
axr.YAxis(2).TickValues=[0.2, 0.4, 0.6];
axr.YAxis(1).Color='b';
axr.YAxis(2).Color='r';
for i=1:length(axr.XTickLabel)
    axr.XTickLabel{i}=['$' axr.XTickLabel{i} '$'];
    axr.XTick(i)=axr.XTick(i);
end
for i=1:length(axr.YTickLabel)
    axr.YTickLabel{i}=['$' axr.YTickLabel{i} '$'];
    axr.YTick(i)=axr.YTick(i);
end
ylabel('$\angle\widetilde{T}$ (rad)', 'Interpreter', 'latex');

title('$z=10$ mm, $y=0$ mm');
xlabel('x (mm)');
leg=legend([h1, h2, h3, h4, h5, h6],...
    'Src. 1 $x,y$', 'Src. 1 $x$',...
    'Det. A $x,y$', 'Det. A $x$',...
    'Det. B $x,y$', 'Det. B $x$',...
    'Location', 'northeast',...
    'NumColumns', 1);
set(leg, 'Position', leg.Position+[0, 0.095, 0, 0]);

h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',5);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 12.5]]);
figName=['exampleT'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');