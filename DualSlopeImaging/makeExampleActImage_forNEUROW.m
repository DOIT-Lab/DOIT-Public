%% Setup
clear; home;

filenames={...
    '211221/HEX1p4_vis_211221';...
    '211222_1/HEX1p4_vis_211222';...
    '211222_2/HEX1p4_vis_211222_2';...
    };

FileInd=1;

ax=[];
filename=filenames{FileInd};
load([filename '_analOutputB.mat']);

Iz_armt=[0, -36.14];

%% Plot O-D
dataTyps={'SDI', 'SDP', 'DSI', 'DSP'};
% dataTyps={'DSI', 'DSP'};
% tits={'Dual-Slope Intensity', 'Dual-Slope Phase'};
% tits={'Single-Distance Intensity', 'Single-Distance Phase',...
%     'Dual-Slope Intensity', 'Dual-Slope Phase'};
tits={'SDI', 'SDP',...
    'DSI', 'DSP'};
useDT=false(size(dataTyps));
for DTind=1:length(dataTyps)
    MTnm=sprintf('%s', dataTyps{DTind}(1:2));
    DTnm=sprintf('%s', dataTyps{DTind}(3));
    
    eval(sprintf(...
        'r_xy=sen.%s.rxy_%s;',...
        MTnm, DTnm));
        
    if ~isempty(r_xy)
        useDT(DTind)=true;
    end
end
sp=[1, sum(useDT)];

col=jet;
col(1, :)=[0, 0, 0];
col(end, :)=[1, 1, 1];
cl=[0, 2.25];

layInd=2;
figure(10000); clf;
h=[];
n=0;
clear ax;
for DTind=1:length(dataTyps)
    MTnm=sprintf('%s', dataTyps{DTind}(1:2));
    DTnm=sprintf('%s', dataTyps{DTind}(3));

    eval(sprintf(...
        'r_xy=sen.%s.rxy_%s;',...
        MTnm, DTnm));

    if ~isempty(r_xy)    
        n=n+1;
        dx=min(diff(unique(r_xy(:, 1))));
        dy=min(diff(unique(r_xy(:, 2))));
        [XX, YY]=meshgrid(min(r_xy(:, 1)):dx:max(r_xy(:, 1)),...
            min(r_xy(:, 2)):dy:max(r_xy(:, 2)));


        hTMP=recon.(MTnm).(DTnm).hOO(:, layInd).*...
            recon.(MTnm).(DTnm).hDD(:, layInd);
        RmAtmp=recon.(MTnm).(DTnm).RmA_OO(:, layInd)-...
            recon.(MTnm).(DTnm).RmA_DD(:, layInd);

        ZZ=griddata(r_xy(:, 1), r_xy(:, 2),...
            -hTMP.*RmAtmp,...
            XX, YY);

        subaxis(sp(1), sp(2), n, 's', 0.01,...
            'ml', 0.08, 'mr', 0.13, 'mb', 0.0); colormap(col);
        pcolor(XX, YY, ZZ); shading flat; hold on;
        plot(Iz_armt(1), Iz_armt(2), 'w.');
%         text(Iz_armt(1)+1, Iz_armt(2), 'Approx. Iz',...
%             'color', 'w');
        hold off;
        caxis(cl);
        ylabel('y (mm)');
        xlabel('x (mm)');
        axis equal tight;
        
        ax(n)=gca;
        if n>1
            ax(n).YTickLabel={};
            ylabel('');
        end
        if n==length(dataTyps)
            cb=colorbar('east');
            cb.Position=cb.Position+[0.04, 0, 0, 0];
            cb.YAxisLocation='right';
            ylabel(cb, sprintf(['\\Delta%s_{Stimulus}-\\Delta%s_{Rest} (\\muM)\n'...
                'where \\Delta%s=\\DeltaO-\\DeltaD'],...
                char(208), char(208), char(208)));
        end

        title(tits{DTind});
    end
end
linkaxes(ax, 'xy');

% sgtitle(sprintf(['Map of Hemoglobin Difference Increase During Stimulation'...
%     '\nMasked by Pixels with Significant Activation']));

h=gcf;
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1);
set(findall(gcf,'-property','Markersize'),'Markersize',3);
set(gcf, 'Units', 'centimeters', 'OuterPosition', [2, 3, [16, 8]]);
figName=['exampleActImage_NEUROW'];
drawnow;
pause(1);
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');