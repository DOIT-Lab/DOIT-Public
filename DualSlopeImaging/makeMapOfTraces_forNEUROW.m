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

DSinds_7x7sp=NaN(7*7, 1);
DSinds_7x7sp(4)=6; % Row1
DSinds_7x7sp(10:12)=[8, 28, 7]; % Row2
DSinds_7x7sp(15:21)=[13, 9, 30, 29, 27, 2, 1]; % Row3
DSinds_7x7sp([22:24, 26:28])=[12, 11, 10, 3, 4, 5]; % Row4
DSinds_7x7sp([29, 31:33, 35])=[15, 24, 25, 26, 23]; % Row5
DSinds_7x7sp([36, 39, 42])=[14, 18, 22]; % Row6
DSinds_7x7sp(44:48)=[16, 17, 19, 20, 21]; % Row7

% for 
FileInd=1; ...:length(filenames)
ax=[];
figure(1000*FileInd); clf;
filename=filenames{FileInd};
load([filename '_analOutputA.mat']);

DScols=hsv(size(DS.rhos, 2));

% Delete Bad Sets
if FileInd==1
    DS.dO_P(:, 16)=NaN(size(DS.dO_P(:, 16)));
end
if FileInd==2
    DS.dO_P(:, 13)=NaN(size(DS.dO_P(:, 13)));
end

for DSind=1:size(DS.rhos, 2)
    
    dtAct=15; %sec
    dtPostAct=30; %sec
    t0Act=t(diff(act)==1);
    
    % fLP=fHPnoise; %Hz
    fLP=0.1; %Hz
    
    optsFold.DTbools=[true, false];
    optsFold.nDT=10;
    
    Cz_armt=mean(armt.rSrc(4:5, 1:2));
    AllSrc_Cz=armt.rSrc(:, 1:2)-Cz_armt;
    AllDet_Cz=armt.rDet(:, 1:2)-Cz_armt;

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
            
            % Delete Bad Sets
            if FileInd==1
                if any([DSind==21, DSind==13, DSind==28, DSind==27,...
                        DSind==7, DSind==1, DSind==16])
                    SD.dO_P(:, SDind)=NaN(size(SD.dO_P(:, SDind)));
                end
            end
            
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


            SDlin_x(n, :)=[AllSrc_Cz(Sind, 1), AllDet_Cz(Dind, 1)];
            SDlin_y(n, :)=[AllSrc_Cz(Sind, 2), AllDet_Cz(Dind, 2)];
        end
    end

    %% Plot Folded Blood
    spInd=find(DSinds_7x7sp==DSind);

    figure(1000*FileInd);
    subaxis(7, 7, spInd, 'sh', 0.01, 'sv', 0.01,...
        'mt', 0.05, 'mb', 0.01, 'mr', 0.01, 'ml', 0.01);
    % SD
    OdataSD=Ofold_SDI;
    DdataSD=Dfold_SDI;
    hold on;
    patch([0, 1, 1, 0]*dtAct, [-1, -1, 1.6, 1.6], DScols(DSind, :),...
        'EdgeColor', 'none', 'FaceAlpha', 0.25);
    for j=1:size(OdataSD, 2)
        [~, k]=min(abs(rhos(j)-[25, 35]));

        if k==1 % Skip short
            continue;
        end

        plot(tFold, OdataSD(:, j), '--',...
            'color', IOcol);
        plot(tFold, DdataSD(:, j), '--',...
            'color', IDcol);
    end
    OdataSD=Ofold_SDP;
    DdataSD=Dfold_SDP;
    for j=1:size(OdataSD, 2)
        [~, k]=min(abs(rhos(j)-[25, 35]));

        if k==1 % Skip short
            continue;
        end

        plot(tFold, OdataSD(:, j), '--',...
            'color', POcol);
        plot(tFold, DdataSD(:, j), '--',...
            'color', PDcol);
    end
    % DS
    Odata=Ofold_DSI;
    Ddata=Dfold_DSI;
    plot(tFold, Odata, '-', 'color', IOcol);
    plot(tFold, Ddata, '-', 'color', IDcol);
    Odata=Ofold_DSP;
    Ddata=Dfold_DSP;
    plot(tFold, Odata, '-', 'color', POcol);
    plot(tFold, Ddata, '-', 'color', PDcol);
    
    yl=ylim;
    plot([1, 1]*dtAct, [-1e6, 1e6], '-.', 'color', [0.5, 0.5, 0.5]);
%     plot([1, 1]*dtAct, [-1e6, 1e6], '-.', 'color', DScols(DSind, :));
    ylim(yl);
    xlim([tFold(1), tFold(end)]);
    hold off;
    axis off;
    
    if spInd==4
        title('Pz');
    elseif spInd==46
        title('Iz');
    end
    
    ax(DSind)=gca;
end

%% Finish Fig
subaxis(7, 7, 25, 'sh', 0.01, 'sv', 0.01);
plot([0, 10]+17, [0, 0], '-k'); hold on;
plot([0, 0]+17, [0, 1], '-k');
text(0+17, -0.3, '10 sec'); 
text(1+17, 0.8, '1 \muM'); hold off;
axis off;
ax(end+1)=gca;

subaxis(7, 7, [1, 2, 3], 'sh', 0.01, 'sv', 0.01);
p1=plot(NaN, NaN, '-k'); hold on;
p2=plot(NaN, NaN, ':', 'color', [0, 0, 0]);
p3=plot(NaN, NaN, '--', 'color', [0, 0, 0]);
p4=plot(NaN, NaN, '-.', 'color', [0, 0, 0]);
leg=legend([p1, p2, p3, p4],...
    'Dual-Slope (25 mm & 37 mm)', 'Single-Distance (25 mm)', 'Single-Distance (37 mm)',...
    'Stimulus End', 'location', 'east');
leg.Position=leg.Position+[-0.005, 0.005, 0, 0];
axis off;

subaxis(7, 7, [5, 6, 7], 'sh', 0.01, 'sv', 0.01);
p5=fill(NaN, NaN, IOcol); hold on;
p6=fill(NaN, NaN, IDcol);
p7=fill(NaN, NaN, POcol);
p8=fill(NaN, NaN, PDcol);
legend([p5, p6, p7, p8],...
    'Intensity \DeltaO', 'Intensity \DeltaD',...
    'Phase \DeltaO', 'Phase \DeltaD');
hold off;
axis off;

linkaxes(ax, 'xy');

%     sgtitle('Dual-Slope Folded Average');
%%
h=gcf;
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1);
set(findall(gcf,'-property','Markersize'),'Markersize',3);
set(gcf, 'Units', 'centimeters', 'OuterPosition', [25, 35, [16, 16]]);
figName=sprintf('traceMap%d_forNEUROW', FileInd);
drawnow;
pause(1);
exportgraphics(h, [figName '.png']);
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');
% end