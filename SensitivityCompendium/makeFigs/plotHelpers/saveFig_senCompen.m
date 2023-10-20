% EXAMPLE PARAMS
% svFig_params.doLatex=true;
% svFig_params.doLatexAxes=[true, true, true]; % X, Y, and Z OR theta and r
% svFig_params.Fontsize=8;
% svFig_params.Linewidth=1.5;
% svFig_params.Markersize=3;
% svFig_params.FigSz_cm=[12, 8];
% svFig_params.doExport=true;

if ~exist('svFig', 'var')
    svFig=true;
end

if svFig
    %% Make Figure Name
    if isempty(h.Name)
        figName=sprintf('%d', h.Number);
    else
        figName=sprintf('%s_%d', h.Name, h.Number);
    end
    
    fprintf('Figure\t%s\n', figName);

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
        svFig_params.doExport=true;
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
        if ~isfield(svFig_params, 'doExport')
            svFig_params.doExport=true;
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
        for iterVar1=1:length(axs)
            dirInd=0;
            for dir=['X', 'Y', 'Z']
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    try
                        axs(iterVar1).([dir 'Axis']).Exponent=...
                            axs(iterVar1).([dir 'Axis']).Exponent;
                    catch
                    end
                end
            end
        end
        for iterVar1=1:length(axs)
            dirInd=0;
            for dir=['X', 'Y', 'Z']
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    try
                        axs(iterVar1).([dir 'Axis']).Exponent=...
                            axs(iterVar1).([dir 'Axis']).Exponent;
                        for iterVar2=1:length(axs(iterVar1).([dir 'Axis']))
                            for iterVar3=1:length(axs(iterVar1).([dir 'Axis'])(iterVar2).TickLabels)
                                axs(iterVar1).([dir 'Axis'])(iterVar2).TickValues(iterVar3)=...
                                    axs(iterVar1).([dir 'Axis'])(iterVar2).TickValues(iterVar3);
                                axs(iterVar1).([dir 'Axis'])(iterVar2).TickLabels{iterVar3}=...
                                    ['$'...
                                    axs(iterVar1).([dir 'Axis'])(iterVar2).TickLabels{iterVar3} '$'];
                            end
                        end
                    catch
                    end
                end
            end
        end

        for iterVar1=1:length(axsPol)
            dirInd=0;
            for dir=["R", "Theta"]
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    try
                        axsPol(iterVar1).([dir 'Axis']).Exponent=...
                            axsPol(iterVar1).([dir 'Axis']).Exponent;
                    catch
                    end
                end
            end
        end
        for iterVar1=1:length(axsPol)
            dirInd=0;
            for dirStr=["R", "Theta"]
                dir=char(dirStr);
                dirInd=dirInd+1;
                if svFig_params.doLatexAxes(dirInd)
                    axsPol(iterVar1).([dir 'Axis']).Exponent=...
                        axsPol(iterVar1).([dir 'Axis']).Exponent;
                    for iterVar2=1:length(axsPol(iterVar1).([dir 'Axis']))
                        for iterVar3=1:length(axsPol(iterVar1).([dir 'Axis'])(iterVar2).TickLabels)
                            axsPol(iterVar1).([dir 'Axis'])(iterVar2).TickValues(iterVar3)=...
                                axsPol(iterVar1).([dir 'Axis'])(iterVar2).TickValues(iterVar3);
                            axsPol(iterVar1).([dir 'Axis'])(iterVar2).TickLabels{iterVar3}=...
                                ['$'...
                                axsPol(iterVar1).([dir 'Axis'])(iterVar2).TickLabels{iterVar3} '$'];
                        end
                    end
                end
            end
        end
        set(findall(gcf, '-property', 'TickLabelInterpreter'),...
            'TickLabelInterpreter','latex');
    end
    
    if ~isempty(svFig_params.Fontsize)
        set(findall(gcf, '-property', 'FontSize'),...
            'FontSize', svFig_params.Fontsize);
    end
    if ~isempty(svFig_params.Linewidth)
        set(findall(gcf, '-property', 'Linewidth'),...
            'Linewidth', svFig_params.Linewidth);
    end
    if ~isempty(svFig_params.Markersize)
        set(findall(gcf, '-property', 'Markersize'),...
            'Markersize', svFig_params.Markersize);
    end
    
    %% Position Figure
    switch deblank(evalc('!hostname'))
        case 'DOIT-Maxwell'
            clNum=mod(svFig_num, 12);
            rwNum=ceil(svFig_num/12+1e-6)-1;
            figPos=[3+clNum*3, 30+rwNum*3];
        case 'DOIT-Odin'
            clNum=mod(svFig_num, 12);
            rwNum=ceil(svFig_num/12+1e-6)-1;
            figPos=[30+clNum*3, 38+rwNum*3];
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
    
    %% Save Figure
    drawnow;
    fprintf('\tDraw done\n')
    
    if svFig_params.doExport
        id1='MATLAB:print:ContentTypeImageSuggested';
        id2='MATLAB:print:ExportExcludesUI';
        warning('off', id1);
        warning('off', id2);
        exportgraphics(h, ['Figs/' figName '.pdf'], 'ContentType', 'vector');
        warning('on', id1);
        warning('on', id2);
        fprintf('\tExport done\n');
    end
end