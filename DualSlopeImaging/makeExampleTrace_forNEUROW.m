%% Setup
clear; home;

filenames={...
    '211221/HEX1p4_vis_211221';...
    '211222_1/HEX1p4_vis_211222';...
    '211222_2/HEX1p4_vis_211222_2';...
    };

IOcol=[1, 0, 0];
IDcol=[0, 0, 1];
POcol=[0.75, 0, 0.5];
PDcol=[0, 0.5, 0.75];

FileInd=1;
DSinds_all=[30, 18];

yl=[-1.5, 2];

ax=[];
filename=filenames{FileInd};
load([filename '_analOutputA.mat']);

DScols=hsv(size(DS.rhos, 2));

dtAct=15; %sec
dtPostAct=30; %sec
t0Act=t(diff(act)==1);

% fLP=fHPnoise; %Hz
fLP=0.5; %Hz

optsFold.DTbools=[true, false];
optsFold.nDT=10;

Iz_armt=[0, -36.14];
AllSrc_Iz=armt.rSrc(:, 1:2)-Iz_armt;
AllDet_Iz=armt.rDet(:, 1:2)-Iz_armt;

for figInd=1:length(DSinds_all)
    DSind=DSinds_all(figInd);

    %% Do Fold
    %DSI
    i=DSind;
    try
        OsigLP=lowpass(DS.dO_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(DS.dD_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        [tFold, Ofold_DSI, ~, OOfold_DSI]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, Dfold_DSI, ~, DDfold_DSI]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    catch
        Ofold_DSI=NaN(size(tFold));
        Dfold_DSI=NaN(size(tFold));
    end
    
    %DSP
    i=DSind;
    try
        OsigLP=lowpass(DS.dO_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(DS.dD_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        [~, Ofold_DSP, ~, OOfold_DSP]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, Dfold_DSP, ~, DDfold_DSP]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    catch
        Ofold_DSP=NaN(size(tFold));
        Dfold_DSP=NaN(size(tFold));
    end
    
    % Find SD sets for DS set
    SSprsTmp=armt.DSprs(:, :, :, DSind);
    n=0;
    for SSind=1:size(SSprsTmp, 3)
        SDprsTmp=SSprsTmp(:, :, SSind);
        for SDindTMP=1:size(SDprsTmp, 1)
            Sind=SDprsTmp(SDindTMP, 1);
            Dind=SDprsTmp(SDindTMP, 2);
    
            SDind=find(and(armt.SDprs(:, 1)==Sind, armt.SDprs(:, 2)==Dind));
    
            n=n+1;
    
            rhos(n)=SD.rho(SDind);
    
            %SDI
            i=SDind;
            try
                OsigLP=lowpass(SD.dO_I(:, i), fLP, 1/median(diff(t)),...
                    'ImpulseResponse', 'iir');
                DsigLP=lowpass(SD.dD_I(:, i), fLP, 1/median(diff(t)),...
                    'ImpulseResponse', 'iir');
                [~, Ofold_SDI(:, n), ~, OOfold_SDI(:, :, n)]=...
                    foldingAvgASync(t, OsigLP,...
                    t0Act, 0, dtAct+dtPostAct, optsFold);
                [~, Dfold_SDI(:, n), ~, DDfold_SDI(:, :, n)]=...
                    foldingAvgASync(t, DsigLP,...
                    t0Act, 0, dtAct+dtPostAct, optsFold);
            catch
                Ofold_SDI(:, n)=NaN(size(tFold));
                Dfold_SDI(:, n)=NaN(size(tFold));
            end
    
            %SDP
            i=SDind;
            try
                OsigLP=lowpass(SD.dO_P(:, i), fLP, 1/median(diff(t)),...
                    'ImpulseResponse', 'iir');
                DsigLP=lowpass(SD.dD_P(:, i), fLP, 1/median(diff(t)),...
                    'ImpulseResponse', 'iir');
                [~, Ofold_SDP(:, n), ~, OOfold_SDP(:, :, n)]=...
                    foldingAvgASync(t, OsigLP,...
                    t0Act, 0, dtAct+dtPostAct, optsFold);
                [~, Dfold_SDP(:, n), ~, DDfold_SDP(:, :, n)]=...
                    foldingAvgASync(t, DsigLP,...
                    t0Act, 0, dtAct+dtPostAct, optsFold);
            catch
                Ofold_SDP(:, n)=NaN(size(tFold));
                Dfold_SDP(:, n)=NaN(size(tFold));
            end
    
    
            SDlin_x(n, :)=[AllSrc_Iz(Sind, 1), AllDet_Iz(Dind, 1)];
            SDlin_y(n, :)=[AllSrc_Iz(Sind, 2), AllDet_Iz(Dind, 2)];
        end
    end
    
    %% Plot Map
%     figure(1001+figInd-1); clf;
%     hSrc=plot(AllSrc_Iz(:, 1), AllSrc_Iz(:, 2), 'sr'); hold on;
%     hDet=plot(AllDet_Iz(:, 1), AllDet_Iz(:, 2), 'ob'); 
%     for i=1:size(SDlin_x, 1)
%         hEx=plot(SDlin_x(i, :), SDlin_y(i, :), '-g');
%     end
%     plot(0, 0, '.k'); 
%     text(1, 0, 'Approx. Iz');
%     hold off;
%     axis equal;
%     xlim([-65, 65]);
%     ylim([-20, 95]);
%     xlabel('x (mm)');
%     ylabel('y (mm)');
%     legend([hSrc, hDet, hEx], 'Source', 'Detector', 'DS set',...
%         'location', 'northeast', 'NumColumns', 1);
%     title('Dual-Slope Set Location');
%     
%     h=gcf;
%     set(findall(gcf,'-property','FontSize'),'FontSize',18);
%     set(findall(gcf,'-property','Linewidth'),'Linewidth',3);
%     set(findall(gcf,'-property','Markersize'),'Markersize',16);
%     set(gcf, 'Units', 'centimeters', 'OuterPosition', [30, 22, [16, 9]*2]);
%     figName=['exampleSetLoc_' num2str(figInd) '.png'];
%     drawnow;
%     pause(1);
%     exportgraphics(h, figName);
    
    %% Plot Example
    figure(2001+figInd-1); clf;
    
    h1=subaxis(1, 2, 1, 's', 0.01, 'ml', 0.09, 'mb', 0.11, 'mt', 0.04,...
        'mr', 0.01);
    OdataSD=Ofold_SDI;
    DdataSD=Dfold_SDI;
    Odata=Ofold_DSI;
    Ddata=Dfold_DSI;
    hold on;
    patch([0, 1, 1, 0]*dtAct, [-1.5, -1.5, 2, 2], DScols(DSind, :),...
        'EdgeColor', 'none', 'FaceAlpha', 0.25);
    plot(tFold, Odata, '-r');
    plot(tFold, Ddata, '-b');
    for j=1:size(OdataSD, 2)
        [~, k]=min(abs(rhos(j)-[25, 35]));
    
        if k==1
            sym=':';
        elseif k==2
            sym='--';
        end
    
        plot(tFold, OdataSD(:, j), sym,...
            'color', IOcol);
        plot(tFold, DdataSD(:, j), sym,...
            'color', IDcol);
    end
%     yl=ylim;
    plot([1, 1]*dtAct, [-1e6, 1e6], '-.k');
    ylim(yl);
    xlim([tFold(1), tFold(end)]);
    p1=plot(NaN, NaN, '-k');
    p2=plot(NaN, NaN, ':', 'color', [0, 0, 0]);
    p3=plot(NaN, NaN, '--', 'color', [0, 0, 0]);
    p4=plot(NaN, NaN, '-.', 'color', [0, 0, 0]);
    legend([p1, p2, p3, p4],...
        'Dual-Slope (25 mm & 37 mm)', 'Single-Distance (25 mm)', 'Single-Distance (37 mm)',...
        'Stimulus End', 'location', 'northeast');
    hold off;
    xlabel('t (sec)');
    ylabel('\muM');
    % title('Intensity');
    
    
    h2=subaxis(1, 2, 2);
    OdataSD=Ofold_SDP;
    DdataSD=Dfold_SDP;
    Odata=Ofold_DSP;
    Ddata=Dfold_DSP;
    hold on;
    patch([0, 1, 1, 0]*dtAct, [-1.5, -1.5, 2, 2], DScols(DSind, :),...
        'EdgeColor', 'none', 'FaceAlpha', 0.25);
    plot(tFold, Odata, '-', 'color', POcol);
    plot(tFold, Ddata, '-', 'color', PDcol);
    for j=1:size(OdataSD, 2)
        [~, k]=min(abs(rhos(j)-[25, 35]));
    
        if k==1
            sym=':';
        elseif k==2
            sym='--';
        end
    
        plot(tFold, OdataSD(:, j), sym,...
            'color', POcol);
        plot(tFold, DdataSD(:, j), sym,...
            'color', PDcol);
    end
%     yl=ylim;
    plot([1, 1]*dtAct, [-1e6, 1e6], '-.k');
    ylim(yl);
    xlim([tFold(1), tFold(end)]);
    p5=fill(NaN, NaN, IOcol);
    p6=fill(NaN, NaN, IDcol);
    p7=fill(NaN, NaN, POcol);
    p8=fill(NaN, NaN, PDcol);
    leg=legend([p5, p6, p7, p8],...
        'Intensity \DeltaO', 'Intensity \DeltaD',...
        'Phase \DeltaO', 'Phase \DeltaD', 'location', 'northeast');
    leg.Position=leg.Position+[0.01, 0, 0, 0];
    hold off;
    xlabel('t (sec)');
    % ylabel('\muM');
    ax=gca;
    ax.YTickLabel={};
    % title('Phase');
    
    linkaxes([h1, h2], 'xy');
    %%
    h=gcf;
    set(findall(gcf,'-property','FontSize'),'FontSize',10);
    set(findall(gcf,'-property','Linewidth'),'Linewidth',1);
    set(findall(gcf,'-property','Markersize'),'Markersize',3);
    set(gcf, 'Units', 'centimeters', 'OuterPosition', [25, 35, [16, 12]]);
    figName=['exampleSetTrace_NEUROW_' num2str(figInd)];
    drawnow;
    pause(1);
    exportgraphics(h, [figName '.png']);
    exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');
end