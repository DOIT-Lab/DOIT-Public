function [Gamma_grid, Delta_grid]=GammaDelta(SpS, sz_grid, dr)
%armt=findSSsets(armt, arry)
% Giles Blaney Summer 2020
% 
% G. Blaney, A. Sassaroli, and S. Fantini, “Design of a source-detector 
% array for dual-slope diffuse optical imaging,” Review of Scientific 
% Instruments, https://doi.org/10.1063/5.0015512.
% 
%   Inputs:
%       SpS     - n X n array of pinv(S)*S
%       sz_grid - 1 X 3 array with size of grid coordinate system 
%                 [n_y, n_x, n_z]
%       dr      - 1 X 3 array with gird spacing [dy, dx, dz]
%   Outputs:
%       Gamma_grid - sz_grid X 3 array of FWHM resolutions in all axis
%                    directions
%       Delta_grid - sz_grid X 3 array of localization errors in all axis
%                    directions

    initVar=NaN([sz_grid, 3]);
    Gamma_grid=initVar;
    Delta_grid=initVar;
    clear initVar;
    for i=1:3
        yesInd=i+3;
        noInds=4:6;
        noInds(noInds==yesInd)=[];

        switch i
            case 1
                [~, R, ~]=meshgrid(1:sz_grid(1), 1:sz_grid(2), 1:sz_grid(3));
            case 2
                [R, ~, ~]=meshgrid(1:sz_grid(1), 1:sz_grid(2), 1:sz_grid(3));
            case 3
                [~, ~, R]=meshgrid(1:sz_grid(1), 1:sz_grid(2), 1:sz_grid(3));
            otherwise
                R=[];
        end

        Xm_grid=sum(sum(reshape(SpS, [sz_grid, sz_grid]),...
            noInds(1)), noInds(2));

        Gamma_grid(:, :, :, i)=sum(Xm_grid>=(max(Xm_grid, [], yesInd)/2),...
            yesInd)*dr(i);
        
        [~, maxInds]=max(Xm_grid, [], yesInd);
        Delta_grid(:, :, :, i)=(maxInds-R)*dr(i);

    end
end

