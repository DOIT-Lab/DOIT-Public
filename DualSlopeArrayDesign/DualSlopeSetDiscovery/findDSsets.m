function armt=findDSsets(armt, arry)
%armt=findSSsets(armt, arry)
% Giles Blaney Summer 2020
% 
% G. Blaney, A. Sassaroli, and S. Fantini, “Design of a source-detector 
% array for dual-slope diffuse optical imaging,” Review of Scientific 
% Instruments, https://doi.org/10.1063/5.0015512.
% 
%   Inputs:
%       armt - Struct with the following fields:
%              - rSrc:  n_src X 3 array of source coordinates [x, y, z;...]
%              - rDet:  n_det X 3 array of detector coordinates [x, y, z;...]
%              - SSprs: 2 X 2 X n_SSprs array of single slope pairs
%                       SSprs(:, :, n_SSprs)=[sInd1, dInd1; sInd2, dInd2]
%       arry - Struct with the following fields:
%              - lTol: Double defining differnece between SSs distance
%                      difference tolerance
%   Outputs:
%       armt - Same as input armt struct with the following added fields:
%              - DSprs: 2 X 2 X 2 X n_DSprs array of dual slope pairs
%                       SSprs(:, :, 1, n)=[sInd11, dInd11; sInd12, dInd12]
%                       SSprs(:, :, 2, n)=[sInd21, dInd21; sInd22, dInd22]

    armt.DSprs=[];
    sdAllHist={};
    n=0;
    for ss1Ind=1:size(armt.SSprs, 3)
        for ss2Ind=1:size(armt.SSprs, 3)
            if ss1Ind>=ss2Ind
                continue;
            end

            s11Ind=armt.SSprs(1, 1, ss1Ind);
            d11Ind=armt.SSprs(1, 2, ss1Ind);
            s21Ind=armt.SSprs(1, 1, ss2Ind);
            d21Ind=armt.SSprs(1, 2, ss2Ind);
            s12Ind=armt.SSprs(2, 1, ss1Ind);
            d12Ind=armt.SSprs(2, 2, ss1Ind);
            s22Ind=armt.SSprs(2, 1, ss2Ind);
            d22Ind=armt.SSprs(2, 2, ss2Ind);

            if length(unique([s11Ind, s21Ind, s12Ind, s22Ind]))==2 && ...
                length(unique([d11Ind, d21Ind, d12Ind, d22Ind]))==2

                sd11=[s11Ind, d11Ind];
                sd21=[s21Ind, d21Ind];
                sd12=[s12Ind, d12Ind];
                sd22=[s22Ind, d22Ind];

                sdAll=[sd11; sd21; sd12; sd22];

                if sum(size(sdAll)==size(unique(sdAll, 'rows')))==2
                    r11=norm(armt.rSrc(s11Ind, :)-armt.rDet(d11Ind, :));
                    r21=norm(armt.rSrc(s21Ind, :)-armt.rDet(d21Ind, :));
                    r12=norm(armt.rSrc(s12Ind, :)-armt.rDet(d12Ind, :));
                    r22=norm(armt.rSrc(s22Ind, :)-armt.rDet(d22Ind, :));

                    [~, ord1]=sort([r11, r12]); %Order of slope 1
                    [~, ord2]=sort([r21, r22]); %Order of slope 2

                    s1Inds=[s11Ind, s12Ind];
                    s1Inds=s1Inds(ord1); %Ordered slope 1 sources
                    s2Inds=[s21Ind, s22Ind];
                    s2Inds=s2Inds(ord2); %Ordered slope 2 sources

                    d1Inds=[d11Ind, d12Ind];
                    d1Inds=d1Inds(ord1); %Ordered slope 1 detectors
                    d2Inds=[d21Ind, d22Ind];
                    d2Inds=d2Inds(ord2); %Ordered slope 2 detectors

                    if sum([s1Inds==s2Inds, d1Inds==d2Inds])>0
                        continue;
                    end

                    if abs(abs(r12-r11)-abs(r22-r21))<=arry.lTol
                        newSet=true;
                        for i=1:length(sdAllHist)
                            if size(unique([sdAllHist{i}; sdAll], 'rows'))<=4
                                newSet=false;
                                break;
                            end
                        end

                        if newSet
                            n=n+1;
                            sdAllHist{n}=sdAll;

                            armt.DSprs(:, :, 1, n)=armt.SSprs(:, :, ss1Ind);
                            armt.DSprs(:, :, 2, n)=armt.SSprs(:, :, ss2Ind);
                        end
                    end
                end
            end
        end
    end
end