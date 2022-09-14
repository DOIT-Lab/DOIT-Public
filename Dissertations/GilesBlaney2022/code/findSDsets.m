function armt=findSDsets(armt, arry)
%armt=findSDsets(armt, arry)
% Giles Blaney
% Inputs:
%   armt - Struct with the following fields:
%          - rSrc: n_src X 3 array of source coordinates:
%              [x, y, z;...]
%          - rDet: n_det X 3 array of detector coordinates:
%              [x, y, z;...]
%   arry - Struct with the following fields:
%          - rRng: 1 X 2 array defining single-distance distance range:
%              [rho_min, rho_max]
% Outputs:
%   armt - Same as input armt struct with the following added fields:
%          - SDprs: n_SDprs X 2 array of single distance pairs:
%              [sInd1, dInd1;...
%                sInd2, dInd2; ...
%                ... ;sIndn, dIndn]

    armt.SDprs=[];
    r_all=NaN(size(armt.rSrc, 1), size(armt.rDet, 1));

    for sInd=1:size(armt.rSrc, 1)
        for dInd=1:size(armt.rDet, 1)
            % Check distance range
            % requirement (arry.rRng)
            r=norm(armt.rSrc(sInd, :)-armt.rDet(dInd, :));
            r_all(sInd, dInd)=r;

            if r<arry.rRng(1) || r>arry.rRng(2)
                continue;
            end

            armt.SDprs(end+1, :)=[sInd, dInd];
        end
    end
end