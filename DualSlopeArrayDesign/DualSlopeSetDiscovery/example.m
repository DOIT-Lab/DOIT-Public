%% Setup
clear; home;
% Giles Blaney Summer 2020

% G. Blaney, A. Sassaroli, and S. Fantini, “Design of a source-detector 
% array for dual-slope diffuse optical imaging,” Review of Scientific 
% Instruments, Submitted.

load('Coors_HEXarray.mat');

armt.rSrc=AllSrcs;
armt.rDet=AllDets;

arry.rRng=[0, 40+1e-3]; %mm
arry.lRng=[10-1e-3, inf]; %mm
arry.lTol=1; %mm

doPlot=true;

%% Find Sets
armt=findSDsets(armt, arry);
armt=findSSsets(armt, arry);
armt=findDSsets(armt, arry);

%% Plot
if doPlot
    col=hsv(size(armt.DSprs, 4));
    
    figure(1); clf;
    for i=1:size(armt.DSprs, 4)
        SSpr1=armt.DSprs(:, :, 1, i);
        SSpr2=armt.DSprs(:, :, 2, i);
        SDprs=[SSpr1; SSpr2];

        sInds=unique(SDprs(:, 1));
        dInds=unique(SDprs(:, 2));

        for k=1:size(SDprs, 1)
            SDprs=[SSpr1; SSpr2];

            sInd=SDprs(k, 1);
            dInd=SDprs(k, 2);

            plot([armt.rSrc(sInd, 1), armt.rDet(dInd, 1)],...
                [armt.rSrc(sInd, 2), armt.rDet(dInd, 2)],...
                '-', 'color', col(i, :));
            hold on;
        end
    end
    h1=plot(armt.rSrc(:, 1), armt.rSrc(:, 2), 'sr'); hold on;
    h2=plot(armt.rDet(:, 1), armt.rDet(:, 2), 'ob'); 
    hold off;
    xlabel('x (mm)');
    ylabel('y (mm)');
    axis equal tight;
    legend([h1, h2], 'Source', 'Detector', 'location', 'bestoutside');
end