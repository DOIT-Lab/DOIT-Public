[cl, cols, cTicks, cTickLabs, cCount]=makeSatCols(plotData(:),...
    'preset', plotPreset);

cLab=['$\mathcal{S}_{\mathcal{Y}}(\vec{r})=' ...
    '\frac{\partial \mathcal{Y} / \partial \mu_{a,pert} (\vec{r})}' ...
    '{\partial \mathcal{Y} / \partial \mu_{a,pert,homo}}=' ...
    '\frac{\Delta \mu_{a,\mathcal{Y}}}' ...
    '{\Delta \mu_{a,pert} (\vec{r})}$'];

dx=diff(xl); dy=diff(yl); dz=diff(zl);
Wx=(svFig_params.FigSz_cm(1)-3*spcH)/(1+dy/dx);
Sy=Wx*dy/dx;
Hz=Wx*dz/dx;

svFig_params.FigSz_cm(2)=2*spcV+padV+spcT+Hz+Sy;

Pos(1, :)=[spcH+padH, Hz+2*spcV+padV, Wx, Sy];
Pos(2, :)=[Wx+2*spcH+padH, svFig_params.FigSz_cm(2)-(spcT+Sy), ...
    svFig_params.FigSz_cm(1)-(Wx+3*spcH+padH), ...
    svFig_params.FigSz_cm(2)-(Hz+2*spcV+padV+spcT)];
Pos(3, :)=[spcH+padH, spcV+padV, Wx, Hz];
Pos(4, :)=[Wx+2*spcH+padH, spcV+padV, Sy, Hz];

%% x-y
sp(1)=axes;
set(sp(1), 'Units', 'centimeters');
set(sp(1), 'Position', Pos(1, :));
set(sp(1), 'PositionConstraint', 'innerposition');
set(sp(1), 'DataAspectRatioMode', 'auto');
[planeData, pp]=sliceS(plotData, params, 'z', zsl);
imagesc(pp.horzAx, pp.vertAx, planeData); hold on;
clim(cl);
set(sp(1), 'Colormap', cols);
xlim(xl); ylim(yl);
if strcmp(plotPreset, 'none')
    contour(pp.horzAx, pp.vertAx, planeData, cCount, '--',...
        'Color', [1, 1, 1]*0.5);  hold off;
else
    contour(pp.horzAx, pp.vertAx, planeData, cTicks, '--',...
        'Color', [1, 1, 1]*0.5);  hold off;
end

pos=get(gca, 'Position');
cb=colorbar('location', 'northoutside');
set(cb, 'Ticks', cTicks);
set(sp(1), 'Position', pos);
cb.Units='centimeters';
cb.Position=[spcH, Pos(1, 2)+Pos(1, 4)+spcV, ...
    svFig_params.FigSz_cm(1)-2*spcH, 0.5];
cb.Ticks=cTicks+[1e-6, zeros(1, length(cTicks)-2), -1e-6];
cb.TickLabels=cTickLabs;
if strcmp(plotPreset, 'none')
    if abs(cTicks(end-1)-cTicks(end))<0.003
        cb.Ticks(end-1)=[];
        cb.TickLabels(end-1)=[];
    end
    if abs(cTicks(1)-cTicks(2))<0.003
        cb.Ticks(2)=[];
        cb.TickLabels(2)=[];
    end
end
cb.TickLabelInterpreter='latex';
cb.Ruler.TickLabelRotation=30;
cb.LineWidth=1.5;
cb.Label.Interpreter='latex';
cb.Label.String=cLab;
cb.Label.FontSize=12;

set(sp(1), 'XTickLabel', {});
ylabel(pp.vertNm, 'Interpreter', 'latex');
title(sprintf('(a) $z=%.1f$ mm', zsl), 'Interpreter', 'latex');

%% x-y-z
sp(2)=axes;
set(sp(2), 'Units', 'centimeters');
set(sp(2), 'Position', Pos(2, :));
set(sp(2), 'PositionConstraint', 'outerposition');
set(sp(2), 'DataAspectRatioMode', 'auto');
hold on;
[XX, YY, ZZ]=meshgrid(params.x, params.y, params.z);
if cl(1)<0
    isosurface(XX, YY, ZZ, permute(plotData, [2, 1, 3]), negSurf);
end
isosurface(XX, YY, ZZ, permute(plotData, [2, 1, 3]), posSurf); 
for i=1:size(rd, 1)
    plot3([0, 0, zl(1)/5, 0, -zl(1)/5]+rd(i, 1), ...
        [0, 0, 0, 0, 0], ...
        [0, zl(1), zl(1)-zl(1)/5, zl(1), zl(1)-zl(1)/5], ...
        '-b', ...
        'LineWidth', 1.5);
end
for i=1:size(rs, 1)
    plot3([0, 0, zl(1)/5, 0, -zl(1)/5]+rs(i, 1), ...
        [0, 0, 0, 0, 0], ...
        [zl(1), 0, zl(1)/5, 0, zl(1)/5], '-r', ...
        'LineWidth', 1.5);
end

plot3(xl, [1, 1]*ysl, [1, 1]*zsl, ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xsl, yl, [1, 1]*zsl, ':', 'color', [0.5, 0.5, 0.5]);

plot3(xl, [1, 1]*yl(1), [1, 1]*zsl, ':', 'color', [0.5, 0.5, 0.5]);
plot3(xl, [1, 1]*yl(2), [1, 1]*zsl, ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(1), yl, [1, 1]*zsl, ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(2), yl, [1, 1]*zsl, ':', 'color', [0.5, 0.5, 0.5]);

plot3([1, 1]*xsl, [1, 1]*yl(1), [0, zl(2)], ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xsl, [1, 1]*ysl, [0, zl(2)], ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xsl, [1, 1]*yl(2), [0, zl(2)], ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xsl, yl, [0, 0], ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xsl, yl, [1, 1]*zl(2), ':', 'color', [0.5, 0.5, 0.5]);

plot3(xl, [1, 1]*ysl, [0, 0], ':', 'color', [0.5, 0.5, 0.5]);
plot3(xl, [1, 1]*ysl, [1, 1]*zl(2), ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(1), [1, 1]*ysl, [0, zl(2)], ':', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(2), [1, 1]*ysl, [0, zl(2)], ':', 'color', [0.5, 0.5, 0.5]);

plot3(xl, [1, 1]*yl(1), [0, 0], '-', 'color', [0.5, 0.5, 0.5]);
plot3(xl, [1, 1]*yl(2), [0, 0], '-', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(1), yl, [0, 0], '-', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(2), yl, [0, 0], '-', 'color', [0.5, 0.5, 0.5]);

plot3(xl, [1, 1]*yl(1), [1, 1]*zl(2), '-', 'color', [0.5, 0.5, 0.5]);
plot3(xl, [1, 1]*yl(2), [1, 1]*zl(2), '-', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(1), yl, [1, 1]*zl(2), '-', 'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(2), yl, [1, 1]*zl(2), '-', 'color', [0.5, 0.5, 0.5]);

plot3([1, 1]*xl(1), [1, 1]*yl(1), [0, zl(2)], '-', ...
    'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(2), [1, 1]*yl(1), [0, zl(2)], '-', ...
    'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(1), [1, 1]*yl(2), [0, zl(2)], '-', ...
    'color', [0.5, 0.5, 0.5]);
plot3([1, 1]*xl(2), [1, 1]*yl(2), [0, zl(2)], '-', ...
    'color', [0.5, 0.5, 0.5]);

hold off;
clim(cl);
set(sp(2), 'Colormap', cols);
set(sp(2), 'ZDir', 'reverse');
axis equal tight;
view(25, 15);
camlight;
xlim(xl); ylim(yl); zlim(zl);

xlabel('$x$ (mm)', 'Interpreter', 'latex');
ylabel('$y$ (mm)', 'Interpreter', 'latex');
zlabel('$z$ (mm)', 'Interpreter', 'latex');
title('(b)', 'Interpreter', 'latex');

%% x-z
sp(3)=axes;
set(sp(3), 'Units', 'centimeters');
set(sp(3), 'Position', Pos(3, :));
set(sp(3), 'PositionConstraint', 'innerposition');
set(sp(3), 'DataAspectRatioMode', 'auto');
[planeData, pp]=sliceS(plotData, params, 'y', ysl);
imagesc(pp.horzAx, pp.vertAx, planeData); hold on;
for i=1:size(rd, 1)
    if rd(i, 2)==ysl
        plot([0, 0, zl(1)/5, 0, -zl(1)/5]+rd(i, 1), ...
            [0, zl(1), zl(1)-zl(1)/5, zl(1), zl(1)-zl(1)/5], '-b', ...
            'LineWidth', 1.5);
    end
end
for i=1:size(rs, 1)
    if rs(i, 2)==ysl
        plot([0, 0, zl(1)/5, 0, -zl(1)/5]+rs(i, 1), ...
            [zl(1), 0, zl(1)/5, 0, zl(1)/5], '-r', ...
            'LineWidth', 1.5);
    end
end
yline(0, ':k');
clim(cl);
set(sp(3), 'Colormap', cols);
xlim(xl); ylim(zl);
if strcmp(plotPreset, 'none')
    contour(pp.horzAx, pp.vertAx, planeData, cCount, '--',...
        'Color', [1, 1, 1]*0.5);  hold off;
else
    contour(pp.horzAx, pp.vertAx, planeData, cTicks, '--',...
        'Color', [1, 1, 1]*0.5);  hold off;
end

xlabel(pp.horzNm, 'Interpreter', 'latex');
ylabel(pp.vertNm, 'Interpreter', 'latex');
title(sprintf('(c) $y=%.1f$ mm', ysl), 'Interpreter', 'latex');

%% y-z
sp(4)=axes;
set(sp(4), 'Units', 'centimeters');
set(sp(4), 'Position', Pos(4, :));
set(sp(4), 'PositionConstraint', 'innerposition');
set(sp(4), 'DataAspectRatioMode', 'auto');
[planeData, pp]=sliceS(plotData, params, 'x', xsl);
imagesc(pp.horzAx, pp.vertAx, planeData); hold on;
for i=1:size(rd, 1)
    if rd(i, 1)==xsl
        plot([0, 0, zl(1)/5, 0, -zl(1)/5]+rd(i, 2), ...
            [0, zl(1), zl(1)-zl(1)/5, zl(1), zl(1)-zl(1)/5], '-b', ...
            'LineWidth', 1.5);
    end
end
for i=1:size(rs, 1)
    if rs(i, 1)==xsl
        plot([0, 0, zl(1)/5, 0, -zl(1)/5]+rs(i, 2), ...
            [zl(1), 0, zl(1)/5, 0, zl(1)/5], '-r', ...
            'LineWidth', 1.5);
    end
end
yline(0, ':k');
clim(cl);
set(sp(4), 'Colormap', cols);
xlim(yl); ylim(zl);
if strcmp(plotPreset, 'none')
    contour(pp.horzAx, pp.vertAx, planeData, cCount, '--',...
        'Color', [1, 1, 1]*0.5);  hold off;
else
    contour(pp.horzAx, pp.vertAx, planeData, cTicks, '--',...
        'Color', [1, 1, 1]*0.5);  hold off;
end
xlabel(pp.horzNm, 'Interpreter', 'latex');

set(sp(4), 'YTickLabel', {});
title(sprintf('(d) $x=%.1f$ mm', xsl), 'Interpreter', 'latex');
