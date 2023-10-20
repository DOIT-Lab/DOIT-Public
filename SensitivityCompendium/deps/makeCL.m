function [cl, cols]=makeCL(x, NVA)
% [cl, cols, cTicks]=makeCL(x)
% 
% Giles Blaney Ph.D. Summer 2023
% 
% Make saturated colormap and limits for a given array 
    
    arguments
        x double;

        NVA.qants=[0.05, 0.95];
    end
    
    cols=jet(100);
    cols(1, :)=0;
    cols(end, :)=1;
    
    cl=[quantile(x(:), NVA.qants(1)), quantile(x(:), NVA.qants(2))];
end