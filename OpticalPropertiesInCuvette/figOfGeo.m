%% Setup
clear; home;

rSrc=[5.5, 5, 0];
rDetA=[16.5, 5, 10];
rDetB=[28.5, 5, 10];

L=[45, 10, 10];

figure(10); clf;
cols=lines(7);

[XX, YY, ZZ]=meshgrid(...
    [-1/2, 1/2]*L(2), [0, 1]*L(3), [-1/2, 1/2]*L(1));

%% Plot z=0 mm
subaxis(1, 3, 1, 's', 0.15,...
    'ml', 0.06, 'mr', 0.06);
hold on;
% Src. 1
plot([1, 1]*rSrc(2)-L(2)/2, [-1, 1]+rSrc(1)-L(1)/2,...
    '-k');
plot([-1, 1]+rSrc(2)-L(2)/2, [1, 1]*rSrc(1)-L(1)/2,...
    '-k');
text(rSrc(2)-L(2)/2-12, rSrc(1)-L(1)/2, 'Src. 1');

% Src. 2
plot([1, 1]*rSrc(2)-L(2)/2, -([-1, 1]+rSrc(1)-L(1)/2),...
    '-', 'Color', [0.5, 0.5, 0.5]);
plot([-1, 1]+rSrc(2)-L(2)/2, -([1, 1]*rSrc(1)-L(1)/2),...
    '-', 'Color', [0.5, 0.5, 0.5]);
text(rSrc(2)-L(2)/2-12, -(rSrc(1)-L(1)/2), 'Src. 2',...
    'Color', [0.5, 0.5, 0.5]);

% Box
plot([-5, 5, 5, -5, -5], [-1, -1, 1, 1, -1]*L(1)/2, '-k');

hold off;
ax2=gca;
for i=1:length(ax2.XTickLabel)
    ax2.XTickLabel{i}=['$' ax2.XTickLabel{i} '$'];
end
for i=1:length(ax2.YTickLabel)
    ax2.YTickLabel{i}=['$' ax2.YTickLabel{i} '$'];
    ax2.YTick(i)=ax2.YTick(i);
end
for i=1:length(ax2.ZTickLabel)
    ax2.ZTickLabel{i}=['$' ax2.ZTickLabel{i} '$'];
end

set(gca, 'YAxisLocation', 'right');

xlabel('$y$ (mm)', 'Interpreter', 'latex');
ylabel('$x$ (mm)', 'Interpreter', 'latex');

title(sprintf('$z=0$ mm\n\\textbf{(a)}'),...
    'Interpreter', 'latex');

axis equal;
xlim([min(XX(:)), max(XX(:))]);
ylim([min(ZZ(:)), max(ZZ(:))]);

%% Plot 3D
subaxis(1, 3, 2);
hold on;

% Src. 1
plot3([1, 1]*rSrc(2)-L(2)/2, [1, 1]*rSrc(3), [-1, 1]+rSrc(1)-L(1)/2,...
    '-k');
plot3([-1, 1]+rSrc(2)-L(2)/2, [1, 1]*rSrc(3), [1, 1]*rSrc(1)-L(1)/2,...
    '-k');

% Src. 2
plot3([1, 1]*rSrc(2)-L(2)/2, [1, 1]*rSrc(3), -([-1, 1]+rSrc(1)-L(1)/2),...
    '-', 'Color', [0.5, 0.5, 0.5]);
plot3([-1, 1]+rSrc(2)-L(2)/2, [1, 1]*rSrc(3), -([1, 1]*rSrc(1)-L(1)/2),...
    '-', 'Color', [0.5, 0.5, 0.5]);

% Det. A
thetaCirc=linspace(0, 2*pi, 100);
rCirc=1;
xCirc=rCirc*cos(thetaCirc);
yCirc=rCirc*sin(thetaCirc);
plot3(yCirc+rDetA(2)-L(2)/2, ones(size(xCirc))*rDetA(3),...
    xCirc+rDetA(1)-L(1)/2, '-', 'Color', cols(1, :));

% Det. B
thetaCirc=linspace(0, 2*pi, 100);
rCirc=1;
xCirc=rCirc*cos(thetaCirc);
yCirc=rCirc*sin(thetaCirc);
plot3(yCirc+rDetB(2)-L(2)/2, ones(size(xCirc))*rDetB(3),...
    xCirc+rDetB(1)-L(1)/2, '-', 'Color', cols(2, :));

% Box
plot3([5, 5, -5, -5, 5], [0, 10, 10, 0, 0], -ones(1, 5)*L(1)/2, '-k');
plot3([5, 5, -5, -5, 5], [0, 10, 10, 0, 0], ones(1, 5)*L(1)/2, '-k');
plot3([5, 5], [0, 0], [-L(1)/2, L(1)/2], '-k');
plot3([5, 5], [10, 10], [-L(1)/2, L(1)/2], '-k');
plot3([-5, -5], [10, 10], [-L(1)/2, L(1)/2], '-k');
plot3([-5, -5], [0, 0], [-L(1)/2, L(1)/2], '-k');

hold off;
ax=gca;
for i=1:length(ax.XTickLabel)
    ax.XTickLabel{i}=['$' ax.XTickLabel{i} '$'];
end
for i=1:length(ax.YTickLabel)
    ax.YTickLabel{i}=['$' ax.YTickLabel{i} '$'];
end
for i=1:length(ax.ZTickLabel)
    ax.ZTickLabel{i}=['$' ax.ZTickLabel{i} '$'];
end

xlabel('$y$ (mm)', 'Interpreter', 'latex');
ylabel('$z$ (mm)', 'Interpreter', 'latex');
zlabel(sprintf('     $x$ (mm)'), 'Interpreter', 'latex');

title(sprintf('\\textbf{(b)}'),...
    'Interpreter', 'latex');

axis equal;
xlim([min(XX(:)), max(XX(:))]);
ylim([min(YY(:)), max(YY(:))]);
zlim([min(ZZ(:)), max(ZZ(:))]);
view([2, 1, 1]);
camlight('headlight');

%% Plot z=10 mm
subaxis(1, 3, 3);
hold on;

% Det. A
thetaCirc=linspace(0, 2*pi, 100);
rCirc=1;
xCirc=rCirc*cos(thetaCirc);
yCirc=rCirc*sin(thetaCirc);
plot(yCirc+rDetA(2)-L(2)/2, xCirc+rDetA(1)-L(1)/2,...
    '-', 'Color', cols(1, :));
text(rDetA(2)-L(2)/2-5.5, rDetA(1)-L(1)/2, 'Det. A',...
    'Color', cols(1, :));

% Det. B
thetaCirc=linspace(0, 2*pi, 100);
rCirc=1;
xCirc=rCirc*cos(thetaCirc);
yCirc=rCirc*sin(thetaCirc);
plot(yCirc+rDetB(2)-L(2)/2, xCirc+rDetB(1)-L(1)/2,...
    '-', 'Color', cols(2, :));
text(rDetB(2)-L(2)/2-5.5, rDetB(1)-L(1)/2, 'Det. B',...
    'Color', cols(2, :));

% Box
plot([-5, 5, 5, -5, -5], [-1, -1, 1, 1, -1]*L(1)/2, '-k');

hold off;
ax1=gca;
ax1.XDir='reverse';
j=0;
for i=1:length(ax1.XTickLabel)
    ax1.XTick(i)=ax1.XTick(i);
    ax1.XTickLabel{i}=['$' ax1.XTickLabel{i} '$'];
    ax1.XTickLabelRotation=0;
end
for i=1:length(ax1.YTickLabel)
    ax1.YTick(i)=ax1.YTick(i);
    ax1.YTickLabel{i}=['$' ax1.YTickLabel{i} '$'];
end

for i=1:length(ax.XTickLabel)
    ax.XTickLabel{i}=['$' ax.XTickLabel{i} '$'];
end
for i=1:length(ax.YTickLabel)
    ax.YTickLabel{i}=['$' ax.YTickLabel{i} '$'];
end
for i=1:length(ax.ZTickLabel)
    ax.ZTickLabel{i}=['$' ax.ZTickLabel{i} '$'];
end

xlabel('$y$ (mm)', 'Interpreter', 'latex');
ylabel('$x$ (mm)', 'Interpreter', 'latex');

title(sprintf('$z=10$ mm\n\\textbf{(c)}'),...
    'Interpreter', 'latex');

axis equal;
xlim([min(XX(:)), max(XX(:))]);
ylim([min(ZZ(:)), max(ZZ(:))]);

h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',5);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 10]]);
figName=['geo'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');