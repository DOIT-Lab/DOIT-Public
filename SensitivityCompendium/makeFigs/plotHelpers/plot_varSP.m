[cl, cols, cTicks, cTickLabs, cCount]=makeSatCols(plotData(:));

cLab=['$\mathcal{S}_{\mathcal{Y}}(\vec{r})=' ...
    '\frac{\partial \mathcal{Y} / \partial \mu_{a,pert} (\vec{r})}' ...
    '{\partial \mathcal{Y} / \partial \mu_{a,pert,homo}}=' ...
    '\frac{\Delta \mu_{a,\mathcal{Y}}}' ...
    '{\Delta \mu_{a,pert} (\vec{r})}$'];

if diff(xl)<40
    sg=[ceil(varNum/3), 3];
else
    sg=[ceil(varNum/2), 2];
end

for i=1:varNum
    [planeData, pp]=sliceS(plotData(:, :, :, i), params, 'y', ysl);

    subaxis(sg(1), sg(2), i, saNVA{:});
    imagesc(pp.horzAx, pp.vertAx, planeData); hold on;
    for j=1:size(rd_all{i}, 1)
        if rd_all{i}(j, 2)==ysl
            plot([0, 0, zl(1)/5, 0, -zl(1)/5]+rd_all{i}(j, 1), ...
                [0, zl(1), zl(1)-zl(1)/5, zl(1), zl(1)-zl(1)/5], '-b', ...
                'LineWidth', 1.5);
        end
    end
    for j=1:size(rs_all{i}, 1)
        if rs_all{i}(j, 2)==ysl
            plot([0, 0, zl(1)/5, 0, -zl(1)/5]+rs_all{i}(j, 1), ...
                [zl(1), 0, zl(1)/5, 0, zl(1)/5], '-r', ...
                'LineWidth', 1.5);
        end
    end

    contour(pp.horzAx, pp.vertAx, planeData, cCount, '--', ...
        'Color', [1, 1, 1]*0.5); 
    hold off;
    axis equal tight;
    yline(0, ':k');
    clim(cl);
    set(gca, 'Colormap', cols);
    xlim(xl); ylim(zl);
    
    title(sprintf(['(%s) ' titFormat], char(96+i), var_all(i)), ...
        'Interpreter', 'latex');

    if i==2
        pos=get(gca, 'Position');
        cb=colorbar('location', 'northoutside');
        set(cb, 'Ticks', cTicks);
        set(gca, 'Position', pos);
        cb.Units='centimeters';
        if max(abs(cCount))>=1
            CBy=2.3;
        else
            CBy=2.1;
        end
        cb.Position=[ ...
            1, svFig_params.FigSz_cm(2)-CBy, ...
            svFig_params.FigSz_cm(1)-2, 1/2];
        cb.Ticks=cTicks+[1e-6, zeros(1, length(cTicks)-2), -1e-6];
        cb.TickLabels=cTickLabs;
        if abs(cTicks(end-1)-cTicks(end))<0.003
            cb.Ticks(end-1)=[];
            cb.TickLabels(end-1)=[];
        end
        if abs(cTicks(1)-cTicks(2))<0.003
            cb.Ticks(2)=[];
            cb.TickLabels(2)=[];
        end
        cb.TickLabelInterpreter='latex';
        cb.Ruler.TickLabelRotation=30;
        cb.LineWidth=1.5;
        cb.Label.Interpreter='latex';
        cb.Label.String=cLab;
        cb.Label.FontSize=12;
    end
    
    if mod(i, sg(2))~=1
        set(gca, 'YTickLabel', {});
    else
        ylabel(pp.vertNm, 'Interpreter', 'latex');
    end
    if i<=(prod(sg)-sg(2)) && ~(mod(varNum, sg(2))~=0 && i==(varNum-(sg(2)-1)))
        set(gca, 'XTickLabel', {});
    else
        xlabel(pp.horzNm, 'Interpreter', 'latex');
    end
end