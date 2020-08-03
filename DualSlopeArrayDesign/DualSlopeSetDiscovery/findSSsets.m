function armt=findSSsets(armt, arry)
%armt=findSSsets(armt, arry)
% Giles Blaney Summer 2020
% 
% G. Blaney, A. Sassaroli, and S. Fantini, “Design of a source-detector 
% array for dual-slope diffuse optical imaging,” Review of Scientific 
% Instruments, Submitted.
% 
%   Inputs:
%       armt - Struct with the following fields:
%              - rSrc:  n_src X 3 array of source coordinates [x, y, z;...]
%              - rDet:  n_det X 3 array of detector coordinates [x, y, z;...]
%              - SDprs: n_SDprs X 2 array of single distance pairs
%                       [sInd1, dInd1; sInd2, dInd2; ... ;sIndn, dIndn]
%       arry - Struct with the following fields:
%              - lRng: 1 X 2 array defining SS distance difference range 
%                      [DELTArho_min, DELTArho_max]
%   Outputs:
%       armt - Same as input armt struct with the following added fields:
%              - SSprs: 2 X 2 X n_SSprs array of single slope pairs
%                       SSprs(:, :, n_SSprs)=[sInd1, dInd1; sInd2, dInd2]

    armt.SSprs=[];
    n=0;
    for sd1Ind=1:size(armt.SDprs, 1)
        for sd2Ind=1:size(armt.SDprs, 1)
            if sd1Ind>=sd2Ind
                continue;
            end

            s1Ind=armt.SDprs(sd1Ind, 1);
            d1Ind=armt.SDprs(sd1Ind, 2);
            s2Ind=armt.SDprs(sd2Ind, 1);
            d2Ind=armt.SDprs(sd2Ind, 2);

            if xor(s1Ind==s2Ind, d1Ind==d2Ind)
                r1=norm(armt.rSrc(s1Ind, :)-armt.rDet(d1Ind, :));
                r2=norm(armt.rSrc(s2Ind, :)-armt.rDet(d2Ind, :));

                % Check lever arm requirement (arry.lRng)
                if r1~=r2 && and(abs(r1-r2)>=arry.lRng(1),...
                        abs(r1-r2)<=arry.lRng(2))
                    n=n+1;

                    if r1<r2
                        armt.SSprs(:, :, n)=[...
                            s1Ind, d1Ind;...
                            s2Ind, d2Ind];
                    elseif r2<r1
                        armt.SSprs(:, :, n)=[...
                            s2Ind, d2Ind;...
                            s1Ind, d1Ind];
                    end
                end
            end
        end
    end
end