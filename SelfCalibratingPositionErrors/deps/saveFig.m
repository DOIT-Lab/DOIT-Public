% Giles Blaney Ph.D.

% EXAMPLE PARAMS
% svFig_params.doLatex=true;
% svFig_params.doLatexAxes=[true, true, true]; % X, Y, and Z OR theta and r
% svFig_params.Fontsize=8;
% svFig_params.Linewidth=1.5;
% svFig_params.Markersize=3;
% svFig_params.FigSz_cm=[12, 8];
% svFig_params.doPNG=false;
% svFig_params.doEPS=false;
% svFig_params.doPDF=true;
% svFig_params.doPDFraster=false;

% FOR LATEX BEAMER
% Set doLaTex=true for last call of saveFig
% doLaTex_params.doSets=false;
% doLaTex_params.doNums=false;
% doLaTex_params.doAll=true;

%% Begin
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
    
    %% Make figs folder 
    id='MATLAB:MKDIR:DirectoryExists';
    warning('off',id);
    try
        mkdir('figs');
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
        svFig_params.doPNG=false;
        svFig_params.doEPS=false;
        svFig_params.doPDF=true;
        svFig_params.doPDFraster=false;
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
        if ~isfield(svFig_params, 'doPNG')
            svFig_params.doPNG=false;
        end
        if ~isfield(svFig_params, 'doEPS')
            svFig_params.doEPS=false;
        end
        if ~isfield(svFig_params, 'doPDF')
            svFig_params.doPDF=true;
        end
        if ~isfield(svFig_params, 'doPDFraster')
            svFig_params.doPDFraster=false;
        end
    end
    
    %% Set Params
    h=gcf;
    if isempty(h.Name)
        tmp=dbstack;
        h.Name=tmp(2).name;
    end
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
                    for axInd=1:length(axs(iterVar1).([dir 'Axis']))
                        axs(iterVar1).([dir 'Axis'])(axInd).Exponent=0;
                    end
                    % try
                    %     axs(iterVar1).([dir 'Axis']).Exponent=...
                    %         axs(iterVar1).([dir 'Axis']).Exponent;
                    % catch
                    % end
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
        case 'doit-odin'
            clNum=mod(svFig_num, 12);
            rwNum=ceil(svFig_num/12+1e-6)-1;
            figPos=[30+clNum*3, 30+rwNum*3];
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

    % if ~exist('figNames_all', 'var')
    %     figNames_all=string(figName);
    % else
    %     figNames_all=[figNames_all; ...
    %         string(figName)];
    % end
    
    %% Save Figure
    drawnow;
    if svFig_params.doPNG
        exportgraphics(h, ['figs/' figName '.png'], 'ContentType', 'image',...
            'Resolution', 300);
    end
    if svFig_params.doEPS
        exportgraphics(h, ['figs/' figName '.eps'], 'ContentType', 'vector');
    end
    if svFig_params.doPDF
        exportgraphics(h, ['figs/' figName '.pdf'], 'ContentType', 'vector');
    end
    if svFig_params.doPDFraster
        exportgraphics(h, ['figs/' figName '.pdf'], 'ContentType', 'image',...
            'Resolution', 300);
    end
    fprintf('Figure\t%s Saved\n', figName);

    %% Cleanup
    clear dir;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Make LaTex
if ~exist('doLaTex')
    doLaTex=false;
end
if doLaTex
    if ~svFig_params.doPDF && ~svFig_params.doPDFraster
        error('doLaTex requires pdf figures (svFig_params.doPDF=true)')
    end
    
    if ~exist('doLaTex_params', 'var')
        doLaTex_params.doSets=false;
        doLaTex_params.doNums=false;
        doLaTex_params.doAll=true;
    else
        if ~isfield(doLaTex_params, 'doSets')
            doLaTex_params.doSets=false;
        end
        if ~isfield(doLaTex_params, 'doNums')
            doLaTex_params.doNums=false;
        end
        if ~isfield(doLaTex_params, 'doAll')
            doLaTex_params.doAll=true;
        end
    end

    %% Get File List
    pdfFiles=dir('figs/*.pdf');
    
    figNums_all=NaN(size(pdfFiles));
    figSets_all=string(NaN(size(pdfFiles)));
    for i=1:length(pdfFiles)
        iSplit=find(pdfFiles(i).name=='_', 1, 'last');
        if isempty(iSplit)
            iSplit=0;
        end
        
        tmp=split(pdfFiles(i).name((iSplit+1):end), '.');
        figNums_all(i)=str2double(tmp(end-1));
        
        tmp=pdfFiles(i).name(1:(iSplit-1));
        figSets_all(i)=tmp;
    end
    figNums_all=sort(unique(figNums_all));
    figSets_all=unique(figSets_all);
    
    %% LaTex of all Files
    if doLaTex_params.doNums
        id='MATLAB:MKDIR:DirectoryExists';
        warning('off',id);
        try
            mkdir('figs/latex/byNum');
        catch
        end
        warning('on',id);
        
        for j=1:length(figNums_all)
            figNamesNums_all=[];
            for i=1:length(figSets_all)
                if figSets_all(i)==""
                    tmp=string(figNums_all(j));
                else
                    tmp=figSets_all(i)+"_"+figNums_all(j);
                end
                if exist("figs/"+tmp+".pdf", 'file')
                    figNamesNums_all=[figNamesNums_all; ...
                        tmp];
                end
            end
            makeANDcompileLaTex('latex/byNum', ...
                "fig"+figNums_all(j), ...
                figNamesNums_all);
        end
    end
    
    if doLaTex_params.doSets
        id='MATLAB:MKDIR:DirectoryExists';
        warning('off',id);
        try
            mkdir('figs/latex/bySet');
        catch
        end
        warning('on',id);
        
        for i=1:length(figSets_all)
            figNamesSets_all=[];
            for j=1:length(figNums_all)
                if figSets_all(i)==""
                    tmp=string(figNums_all(j));
                else
                    tmp=figSets_all(i)+"_"+figNums_all(j);
                end
                if exist("figs/"+tmp+".pdf", 'file')
                    figNamesSets_all=[figNamesSets_all; ...
                        tmp];
                end
            end
            makeANDcompileLaTex('latex/bySet', ...
                figSets_all(i), ...
                figNamesSets_all);
        end
    end
    
    if doLaTex_params.doAll
        id='MATLAB:MKDIR:DirectoryExists';
        warning('off',id);
        try
            mkdir('figs/latex/allFigs');
        catch
        end
        warning('on',id);

        figNames_all=[];
        for i=1:length(figSets_all)
            for j=1:length(figNums_all)
                if figSets_all(i)==""
                    tmp=string(figNums_all(j));
                else
                    tmp=figSets_all(i)+"_"+figNums_all(j);
                end
                if exist("figs/"+tmp+".pdf", 'file')
                    figNames_all=[figNames_all; ...
                        tmp];
                end
            end
        end
        makeANDcompileLaTex('latex/allFigs', 'allFigs', figNames_all);
    end
    
end
doLaTex=false;

%% Functions
function [] = makeANDcompileLaTex(LaTexFolder, LaTexName, figNames_all)
    
    tit1=deblank(evalc('!hostname'));
    
    wd=strrep(pwd, '\', '/');
    tit2=strrep(strrep(wd, '_', '\_'), '/', ' / ');

    tit3=strrep(LaTexName, '_', '\_');
    
    fid=fopen("figs/"+LaTexFolder+"/"+LaTexName+".tex", 'w+');
    
    fprintf(fid, '\\documentclass[11pt]{beamer}\n');
    fprintf(fid, '\\beamertemplatenavigationsymbolsempty\n');
    fprintf(fid, '\\usepackage{tikz}\n\n');
    
    fprintf(fid, '\\begin{document}\n\n');
    
    fprintf(fid, '\t\\title{%s\\\\%s}\n\n', tit1, tit2);
    fprintf(fid, '\t\\subtitle{%s}\n\n', tit3);
    
    fprintf(fid, '\t\\begin{frame}[c]\n');
    fprintf(fid, '\t\t\\maketitle\n');
    fprintf(fid, '\t\\end{frame}\n\n');
    
    fprintf(fid, '\t\\foreach \\figNm in {');
    for i=1:length(figNames_all)
        fprintf(fid, '%s', figNames_all(i));
        if i~=length(figNames_all)
            fprintf(fid, ',');
        else
            fprintf(fid, '}{\n');
        end
    end
    fprintf(fid, '\t\t\\begin{frame}[c]\n');
    fprintf(fid, '\t\t\t\\centering\n');
    fprintf(fid, ['\t\t\t\\includegraphics' ...
        '[width=\\textwidth,height=\\textheight,keepaspectratio]']);
    fprintf(fid, '{../../\\figNm .pdf}\n');
    fprintf(fid, '\t\t\\end{frame}\n');
    fprintf(fid, '\t}\n\n');
    
    fprintf(fid, '\\end{document}');
    
    fclose(fid);
    
    cd(['figs/' LaTexFolder]);
    system("pdflatex "+LaTexName+".tex");
    cd('..'); cd('..'); cd('..');
end