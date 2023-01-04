% EXAMPLE PARAMS
% svFig_params.doLatex=true;
% svFig_params.doLatexAxes=[true, true, true]; % X, Y, and Z OR theta and r
% svFig_params.Fontsize=8;
% svFig_params.Linewidth=1.5;
% svFig_params.Markersize=3;
% svFig_params.FigSz_cm=[12, 8];

if ~exist('svFig', 'var')
    svFig=true;
end

if svFig
    %% Find number saved since last clear
    if ~exist('svFig_num', 'var')
        svFig_num=0;
    else
        svFig_num=mod(svFig_num+1, 60);
    end
    
    %% Make Figs folder 
    id='MATLAB:MKDIR:DirectoryExists';
    warning('off',id);
    try
        mkdir('Figs');
    catch
    end
    warning('on',id);
    
    %% Set default params if not already set
    if ~exist('svFig_params', 'var')
        svFig_params.doLatex=true;
        svFig_params.doLatexAxes=true(1, 3);
        svFig_params.Fontsize=8;
        svFig_params.Linewidth=1.5;
        svFig_params.Markersize=3;
        svFig_params.FigSz_cm=[12, 8];
    else
        if ~isfield(svFig_params, 'doLatex')
            svFig_params.doLatex=true;
        end
        if ~isfield(svFig_params, 'doLatexAxes')
            svFig_params.doLatexAxes=[true, true, true];
        elseif length(svFig_params.doLatexAxes)<3
            svFig_params.doLatexAxes=[svFig_params.doLatexAxes,...
                true(1, 3-length(length(svFig_params.doLatexAxes)))];
        end
        if ~isfield(svFig_params, 'Fontsize')
            svFig_params.Fontsize=8;
        end
        if ~isfield(svFig_params, 'Linewidth')
            svFig_params.Linewidth=1.5;
        end
        if ~isfield(svFig_params, 'Markersize')
            svFig_params.Markersize=3;
        end
        if ~isfield(svFig_params, 'FigSz_cm')
            svFig_params.FigSz_cm=[12, 8];
        end
    end
    
    %% Set Params
    h=gcf;
    axs=findall(h, 'type', 'axes');
    axsPol=findall(h, 'type', 'polaraxes');
    if svFig_params.doLatex
        set(findall(gcf, '-property', 'Interpreter'),...
            'Interpreter', 'latex');
    end
    if any(svFig_params.doLatexAxes)
        % If Latex we need to add $$ around the axis tick labels
        for i=1:length(axs)
            dirInd=0;
            for dir=['X', 'Y', 'Z']
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    try
                        axs(i).([dir 'Axis']).Exponent=...
                            axs(i).([dir 'Axis']).Exponent;
                    catch
                    end
                end
            end
        end
        for i=1:length(axs)
            dirInd=0;
            for dir=['X', 'Y', 'Z']
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    try
                        axs(i).([dir 'Axis']).Exponent=...
                            axs(i).([dir 'Axis']).Exponent;
                        for j=1:length(axs(i).([dir 'Axis']))
                            for k=1:length(axs(i).([dir 'Axis'])(j).TickLabels)
                                axs(i).([dir 'Axis'])(j).TickValues(k)=...
                                    axs(i).([dir 'Axis'])(j).TickValues(k);
                                axs(i).([dir 'Axis'])(j).TickLabels{k}=...
                                    ['$'...
                                    axs(i).([dir 'Axis'])(j).TickLabels{k} '$'];
                            end
                        end
                    catch
                    end
                end
            end
        end

        for i=1:length(axsPol)
            dirInd=0;
            for dir=["R", "Theta"]
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    try
                        axsPol(i).([dir 'Axis']).Exponent=...
                            axsPol(i).([dir 'Axis']).Exponent;
                    catch
                    end
                end
            end
        end
        for i=1:length(axsPol)
            dirInd=0;
            for dirStr=["R", "Theta"]
                dir=char(dirStr);
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    axsPol(i).([dir 'Axis']).Exponent=...
                        axsPol(i).([dir 'Axis']).Exponent;
                    for j=1:length(axsPol(i).([dir 'Axis']))
                        for k=1:length(axsPol(i).([dir 'Axis'])(j).TickLabels)
                            axsPol(i).([dir 'Axis'])(j).TickValues(k)=...
                                axsPol(i).([dir 'Axis'])(j).TickValues(k);
                            axsPol(i).([dir 'Axis'])(j).TickLabels{k}=...
                                ['$'...
                                axsPol(i).([dir 'Axis'])(j).TickLabels{k} '$'];
                        end
                    end
                end
            end
        end
        set(findall(gcf, '-property', 'TickLabelInterpreter'),...
            'TickLabelInterpreter','latex');
    end

    set(findall(gcf, '-property', 'FontSize'),...
        'FontSize', svFig_params.Fontsize);
    set(findall(gcf, '-property', 'Linewidth'),...
        'Linewidth', svFig_params.Linewidth);
    set(findall(gcf, '-property', 'Markersize'),...
        'Markersize', svFig_params.Markersize);
    
    %% Position Figure
    switch deblank(evalc('!hostname'))
        case 'DOIT-Maxwell'
            clNum=mod(svFig_num, 12);
            rwNum=ceil(svFig_num/12+1e-6)-1;
            figPos=[3+clNum*3, 30+rwNum*3];
        case 'freyr'
            clNum=mod(svFig_num, 10);
            rwNum=ceil(svFig_num/10+1e-6)-1;
            figPos=[1+clNum*3, 1+rwNum*3];
        otherwise
            warning('Hostname not set in saveFig');
            figPos=[0, 0];
    end
    
    set(gcf, 'Units', 'centimeters', 'Innerposition', [figPos,...
        svFig_params.FigSz_cm]);
    
    %% Make Figure Name
    if isempty(h.Name)
        figName=sprintf('%d', h.Number);
    else
        figName=sprintf('%s_%d', h.Name, h.Number);
    end
    
    %% Save Figure
    drawnow;
    exportgraphics(h, ['Figs/' figName '.png'], 'ContentType', 'image',...
        'Resolution', 300);
    exportgraphics(h, ['Figs/' figName '.eps'], 'ContentType', 'vector');
    fprintf('Figure\t%s Saved\n', figName);
end