%PreProcessing for 2 fiber 4 PMT DiFC

clear; close all;
probe_distance = 3;
%% Loading data,Peaks candidates found for Fiber 1: 33

%Allows for selection of folder with .mat data files
saveDir = uigetdir();

% proccodes library (where processing functions are stored)
%Add proccodes path for your own computer before running
% ***************************************************
proccodes_library ="C:\Users\niedr\OneDrive - Northeastern University\Niedre_Lab\Fernando\";
% ***************************************************
addpath(genpath(proccodes_library));

%Getting file name
if ispc
    slash = '\';
else
    slash = '/';
end
index = find(saveDir == slash,1,'last');
stem = saveDir(index+1:end);
if saveDir(end) ~= slash; saveDir = [saveDir slash]; end
fname = [saveDir stem];

rel_thresh = [5 5];
%-----------------------------------------------------------------------------------------------------------------%
%% Running proccessing code
fprintf('Running 2 Fiber 4 Current PMTs preProc...\n\n')

fprintf('Pre-processing %s...\n',stem)

% Load data from both fibers
 load([fname '_DA.mat'], 'time', 'data', 'params')
data_sum1 = data(:,1) + data(:,2);

params1 = params;

 load([fname '_DB.mat'], 'data', 'params')


data_sum2 = data(:,1) + data(:,2);


params2 = params;

% Sampling frequnecy
fs = 1 ./ (time(2) - time(1));

% Format data from both fibers
data = [data_sum1 data_sum2];
params = [params1(1) params2(1)];
params(1).name = [stem ' Fiber 1'];
params(2).name = [stem ' Fiber 2'];
clear data_sum1 data_sum2 params1 params2

% Pre-proc background subtracts data, calculates the noise/ moving peak
% threshhold, and identifies peak candidates. Processes all data in the
% data array
[data_bs, noise, peaks, thresh_curve,in_dat] = preProc(data, time, params, 'RelativeThresh', rel_thresh,'ProminenceFactor',1);
%decrease prominence factor increases peak candidates

%,'Mode','setThresh', 'HardThresh', []
% Calculating SNR based on the estimated run noise, not control noise
snr_fiber_1 = mean(20*log10(peaks(1).pks./noise(1)));
snr_fiber_2 = mean(20*log10(peaks(2).pks./noise(2)));

pk_per_min_fiber_1 = peaks(1).count/time(end);
pk_per_min_fiber_2 = peaks(2).count/time(end);
%Calculating mean background level of the orginal data
mean_background = mean(data);
%-----------------------------------------------------------------------------------------------------------------%
% REMOVE COINCIDENT PEAKS
%-----------------------------------------------------------------------------------------------------------------%

 [peaks(1), peaks(2), coinc_pk_count] =... 
removeCoincPeaks(peaks(1), peaks(2), time, 'CoincidenceWindow',0.03);
 %coinc_pk_count = 0;
plotPeaks(data_bs, time,peaks, thresh_curve, params);
%-----------------------------------------------------------------------------------------------------------------%
% MATCH PEAKS IN DIRECTIONS
%-----------------------------------------------------------------------------------------------------------------%
disp('Matching peaks in the forward direction')
[fwd_peaks(1), fwd_peaks(2), fwd_speed, fwd_score] =...
    matchDirectionalPeaks(peaks(1), peaks(2), time, probe_distance);
plotPeaks(data_bs, time, fwd_peaks, thresh_curve, params,'Direction', 'fwd');
disp('Matching peaks in the reverse direction')
[rev_peaks(2), rev_peaks(1), rev_speed, rev_score] =...
    matchDirectionalPeaks(peaks(2), peaks(1), time, probe_distance);
plotPeaks(data_bs, time, rev_peaks, thresh_curve, params,'Direction', 'rev');
% Save processed data, including parameters for processing
disp('Saving processed data')

plotAllPeaks2(data_bs, time, peaks, fwd_peaks, rev_peaks, thresh_curve, params)

pname = sprintf('%s_proc_relThresh_%g_%g.mat', fname, rel_thresh(1), rel_thresh(2));
save(pname, 'time', 'data_bs', 'params', 'noise','snr_fiber_1','snr_fiber_2', 'peaks','mean_background', 'thresh_curve', 'in_dat');
%-----------------------------------------------------------------------------------------------------------------%
% Display useful results
%-----------------------------------------------------------------------------------------------------------------%
fprintf('\n\n\n--------------------------------------------------------\n')
try
     fprintf('Detected %4.2g peaks per min in Fiber 1\n', pk_per_min_fiber_1);
    fprintf('Detected %4.2g peaks per min in Fiber 2\n', pk_per_min_fiber_2);
     fprintf('Removed a total of %g coincident peaks\n',coinc_pk_count);
    fprintf('Peaks candidates found for Fiber 1: %g\n', peaks(1).count);
     fprintf('Peaks candidates found for Fiber 2: %g\n',peaks(2).count);
     fprintf('Arterial Direction Matched Peaks: %g\n', fwd_peaks(1).count);
     fprintf('Venous Matched Peaks: %g\n', rev_peaks(1).count);
    fprintf('Calculated noise for fiber 1: %.3f nA\n', noise(1));
    fprintf('Calculated noise for fiber 2: %.3f nA\n', noise(2));
     fprintf('Estimated Fiber 1 SNR: %.3f dB\n', snr_fiber_1)
     fprintf('Estimated Fiber 2 SNR: %.3f dB\n', snr_fiber_2)
%     fprintf('Average relative threshold for Fiber 1: %.3f nA\n', mean(thresh_curve(1)));
%     fprintf('Average relative threshold for Fiber 2: %.3f nA\n', mean(thresh_curve(2)));
    fprintf('Average background level for Fiber 1: %.3f nA\n',mean_background(1));
    fprintf('Average background level for Fiber 2: %.3f nA\n',mean_background(2));
     fprintf('Average peak amplitude for Fiber 1: %.3f nA\n', mean(peaks(1).pks));
     fprintf('Average peak amplitude for Fiber 2: %.3f nA\n', mean(peaks(2).pks));
%     fprintf('Peak amplitude standard deviation for Fiber 1: %.3f nA\n', std(peaks(1).pks));
%     fprintf('Peak amplitude standard deviation for Fiber 2: %.3f nA\n', std(peaks(2).pks));
  %   fprintf('Largest peak amplitude for Fiber 1: %.3f nA\n', max(peaks(1).pks));
   %  fprintf('Smallest peak amplitude for Fiber 1: %.3f nA\n', min(peaks(1).pks));
%     fprintf('Largest peak amplitude for Fiber 2: %.3f nA\n', max(peaks(2).pks));
%     fprintf('Smallest peak amplitude for Fiber 2: %.3f nA\n', min(peaks(2).pks));
catch
     fprintf('Peaks found for Fiber 1: %g\n', peaks(1).count);
     fprintf('Peaks found for Fiber 2: %g\n', peaks(2).count);
    fprintf('Calculated noise for fiber 1: %.3f nA\n', noise(1));
    fprintf('Calculated noise for fiber 2: %.3f nA\n', noise(2));
         fprintf('Estimated Fiber 1 SNR: %.3f dB\n', snr_fiber_1)
         fprintf('Estimated Fiber 2 SNR: %.3f dB\n', snr_fiber_2)
%     fprintf('Average relative threshold for Fiber 1: %.3f nA\n', mean(thresh_curve(1)));
%     fprintf('Average relative threshold for Fiber 2: %.3f nA\n', mean(thresh_curve(2)));
    fprintf('Average background level for Fiber 1: %.3f nA\n',mean_background(1));
    fprintf('Average background level for Fiber 2: %.3f nA\n',mean_background(2));
     fprintf('Average peak amplitude for Fiber 1: %.3f nA\n', mean(peaks(1).pks));
     fprintf('Average peak amplitude for Fiber 2: %.3f nA\n', mean(peaks(2).pks));
%     fprintf('Peak amplitude standard deviation for Fiber 1: %.3f nA\n', std(peaks(1).pks));
%     fprintf('Peak amplitude standard deviation for Fiber 2: %.3f nA\n', std(peaks(2).pks));
%     fprintf('Largest peak amplitude for Fiber 1: %.3f nA\n', max(peaks(1).pks));
%     fprintf('Smallest peak amplitude for Fiber 1: %.3f nA\n', min(peaks(1).pks));
%     fprintf('Largest peak amplitude for Fiber 2: %.3f nA\n', max(peaks(2).pks));
%     fprintf('Smallest peak amplitude for Fiber 2: %.3f nA\n', min(peaks(2).pks));
end