function [cl, cols, cTicks, cTickLabs, cCount]=makeSatCols(x, NVA)
    
    arguments
        x double;
        
        NVA.preset string = 'none';

        NVA.clSat0 (1,2) = [-0.02, 0.05];
        NVA.cTicks0Step (1,1) = 0.01;
        NVA.BlkWhtFrac (1,1) = 0.2;

        NVA.numCols (1,1) = 1e3;
        NVA.colMap string = 'turbo';

        NVA.decPnts (1,1) = 2;
        NVA.useSciNot logical = false;
    end
    
    switch NVA.preset
        case 'vox'
            clSat0=[0, 5e-5];
            cTicks0Step=1e-6;
            BlkWhtFrac=0.2;
            numCols=1e3;
            colMap='parula';
            decPnts=5;
            useSciNot=true;
        case 'none'
            clSat0=NVA.clSat0;
            cTicks0Step=NVA.cTicks0Step;
            BlkWhtFrac=NVA.BlkWhtFrac;
            numCols=NVA.numCols;
            colMap=NVA.colMap;
            decPnts=NVA.decPnts;
            useSciNot=NVA.useSciNot;
        otherwise
            error('Unknown preset')
    end
    
    cTicks0=clSat0(1):cTicks0Step:clSat0(2);
    
    xMin=min([x(:); 0]);
    xMax=max(x(:));
    
    clSat=[max([clSat0(1), xMin]), min([clSat0(2), xMax])];
    if clSat(1)==xMin
        cl(1)=clSat(1);
    else
        cl(1)=clSat0(1)-diff(clSat0)*BlkWhtFrac;
    end
    if clSat(2)==xMax
        cl(2)=clSat(2);
    else
        cl(2)=clSat0(2)+diff(clSat0)*BlkWhtFrac;
    end
    
    cSatAx=linspace(clSat0(1), clSat0(2), numCols);
    eval(sprintf('colsSat=%s(length(cSatAx));', colMap));
    indsL=cSatAx<xMin;
    indsU=cSatAx>xMax;
    cSatAx(or(indsL, indsU))=[];
    colsSat(or(indsL, indsU), :)=[];

    cAx=cl(1):median(diff(cSatAx)):cl(2);
    cols=NaN(length(cAx), 3);
    [~, i]=min(abs(cAx-floor(cSatAx(1)*10^(decPnts+1))/10^(decPnts+1) ));
    cols(i:(i+length(cSatAx)), :)=[colsSat; colsSat(end, :)];
    if i>1
        cols(1:(i-1), :)=0;
    end
    if (i+length(cSatAx))<length(cAx)
        cols((i+length(cSatAx)+1):length(cAx), :)=1;
    end
    
    cTicks=cTicks0;
    cTicks(cTicks<cSatAx(1))=[];
    cTicks(cTicks>cSatAx(end))=[];
    if useSciNot
        cCount=unique(round([xMin, cTicks, xMax], decPnts+1));
        cTicks=unique(round([cl(1), cTicks, cl(2)], decPnts));
    else
        cCount=unique(round([xMin, cTicks, xMax], 12));
        cTicks=unique(round([cl(1), cTicks, cl(2)], 12));
    end
    
    if useSciNot
        cTickLabs=cell(size(cCount));
        for i=1:length(cCount)
            if round(cCount(i), decPnts+1)==0
                cTickLabs{i}=sprintf('%.0f', cCount(i));
            else
                cTickLabs{i}=sprintf('%.0e', cCount(i));
                cTickLabs{i}=['$' strrep(cTickLabs{i}, 'e-0', ...
                    '\times 10^{-') '}$'];
            end
        end
    else
        cTickLabs=cell(size(cCount));
        for i=1:length(cCount)
            if i==length(cCount)
                cTickLabs{i}=eval([ ...
                    'sprintf(''%.' num2str(decPnts+1) 'f'', cCount(i));']);
            elseif i==1
                if cCount(i)==0
                    cTickLabs{i}=sprintf('%.0f', cCount(i));
                else
                    cTickLabs{i}=eval([ ...
                    'sprintf(''%.' num2str(decPnts+1) 'f'', cCount(i));']);
                end
            else
                cTickLabs{i}=eval([ ...
                    'sprintf(''%.' num2str(decPnts) 'f'', cCount(i));']);
            end
        end
    end
end