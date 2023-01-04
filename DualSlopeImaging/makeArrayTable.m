%% Setup
clear; home;

filenames={...
    '211221/HEX1p4_vis_211221';...
    '211222_1/HEX1p4_vis_211222';...
    '211222_2/HEX1p4_vis_211222_2';...
    };

DSinds_7x7sp=NaN(7*7, 1);
DSinds_7x7sp(4)=6; % Row1
DSinds_7x7sp(10:12)=[8, 28, 7]; % Row2
DSinds_7x7sp(15:21)=[13, 9, 30, 29, 27, 2, 1]; % Row3
DSinds_7x7sp([22:24, 26:28])=[12, 11, 10, 3, 4, 5]; % Row4
DSinds_7x7sp([29, 31:33, 35])=[15, 24, 25, 26, 23]; % Row5
DSinds_7x7sp([36, 39, 42])=[14, 18, 22]; % Row6
DSinds_7x7sp(44:48)=[16, 17, 19, 20, 21]; % Row7

FileInd=1;
filename=filenames{FileInd};
load([filename '_analOutputA.mat'], 'DS', 'armt');

DScols=hsv(size(DS.rhos, 2));

%% Colors
fprintf('\n');
for i=1:size(DS.rhos, 2)
    DSprs=armt.DSprs(:, :, :, i);
    SSprs1=DSprs(:, :, 1);
    SSprs2=DSprs(:, :, 2);

    srcInds=unique(sort([SSprs1(:, 1); SSprs2(:, 1)]));
    detInds=unique(sort([SSprs1(:, 2); SSprs2(:, 2)]));
    
    figPanInd=find(i==figPanOrd);
    if figPanInd<=26
        figPanLet=[' (' char(96+figPanInd) ')'];
    else 
        figPanLet=['(a' char(96+figPanInd-26) ')'];
    end
    
    fprintf('\\definecolor{DScol%d}{rgb}{%.1f,%.1f,%.1f} ',...
        i, DScols(i, :));

    fprintf('\n');
end

%% Table
figPanOrd=DSinds_7x7sp(~isnan(DSinds_7x7sp));

Iz2Oz=365; %mm
rIz=[0, -40]; %mm
rOz=rIz+[0, 0.1*Iz2Oz];
rPz=rIz+[0, 0.3*Iz2Oz];

fprintf('\n');
for i=1:size(DS.rhos, 2)
    DSprs=armt.DSprs(:, :, :, i);
    SSprs1=DSprs(:, :, 1);
    SSprs2=DSprs(:, :, 2);

    srcInds=unique(sort([SSprs1(:, 1); SSprs2(:, 1)]));
    detInds=unique(sort([SSprs1(:, 2); SSprs2(:, 2)]));
    
    srcPoss=armt.rSrc(srcInds, 1:2);
    detPoss=armt.rDet(detInds, 1:2);

    xCent=mean([srcPoss(:, 1); detPoss(:, 1)]);
    yCent=mean([srcPoss(:, 2); detPoss(:, 2)]);
    
    rhoIz=sqrt((rIz(1)-xCent)^2+(rIz(2)-yCent)^2);
    rhoOz=sqrt((rOz(1)-xCent)^2+(rOz(2)-yCent)^2);
    rhoPz=sqrt((rPz(1)-xCent)^2+(rPz(2)-yCent)^2);

    figPanInd=find(i==figPanOrd);
    if figPanInd<=26
        figPanLet=[' (' char(96+figPanInd) ')'];
    else 
        figPanLet=['(a' char(96+figPanInd-26) ')'];
    end
    
    fprintf('\\rowcolor{DScol%d!25} ', i);
    fprintf('%2.0d & ', i); %DSind
    fprintf('%2.0f & %2.0f & ', xCent, yCent); % Cent Pos
    fprintf('%2.0f & %2.0f & %2.0f & ', rhoIz, rhoOz, rhoPz);
    fprintf('\\texttt{%2.0d} & \\texttt{%c} & \\texttt{%c} & \\texttt{%2.0d} & ',...
        srcInds(1), char(64+detInds(1)), char(64+detInds(2)), srcInds(2)); %Optodes
    fprintf('%s ', figPanLet);

    fprintf('\\\\\n');
end
