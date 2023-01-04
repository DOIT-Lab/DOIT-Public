% Giles Blaney Ph.D. 2022

close all;
addpath('deps/');

A1_look;
fprintf('A1 done\n');

A2_look_blood;
fprintf('A2 done\n');

A3_look_fold;
fprintf('A3 done\n');

B1_makeSen;
fprintf('B1 done\n');

B2_doRecon;
fprintf('B2 done\n');

B3_ttest;
fprintf('B3 done\n');

rmpath('deps/');