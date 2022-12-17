clc
clear all
close all

%%% Description %%%%
% This script load MDSI inference results and do statistical analysis to
% identify significant causal interactions during oddball task.

Current_dir=pwd;
Model_Dir=fullfile(Current_dir,'/', 'model');

%add relevant paths
idcs = strfind(Current_dir,'/');
utils_dir = Current_dir(1:idcs(end-2)-1);
addpath(fullfile(utils_dir,'/','utils'));

pval=1; %0.01

% Window length: 2sec.
window=2;
win2sec=load(fullfile(Model_Dir, '/', sprintf('Stats_GCaMP_oddball_win%d_pval0p0%d.mat', window, pval)));

% Plot causal interaction of each task conditions
figure;
subplot(1,3,1);
imagesc(win2sec.AC_sig(:,:,1));
ylabel('To');
xlabel('From');
title('Causal interactions during control stimulus');
caxis([-0.15 0.15]);
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
set(gca, 'XTick', 1:4, ...                             % Change the axes tick marks
         'XTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...  %   and tick labels
         'YTick', 1:4, ...
         'YTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...
         'TickLength', [0 0]);
colorbar('southoutside');
colormap(mycolormap);

subplot(1,3,2);
imagesc(win2sec.AC_sig(:,:,2));
ylabel('To');
xlabel('From');
title('Causal interactions during oddball');
caxis([-0.15 0.15]);
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
set(gca, 'XTick', 1:4, ...                             % Change the axes tick marks
         'XTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...  %   and tick labels
         'YTick', 1:4, ...
         'YTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...
         'TickLength', [0 0]);
colorbar('southoutside');
colormap(mycolormap);

subplot(1,3,3);
imagesc(win2sec.sig_t_mtx_ONvsOFF);
ylabel('To');
xlabel('From');
title('Causal interactions that distinguish oddball from control');
caxis([-8 8]);
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
set(gca, 'XTick', 1:4, ...                             % Change the axes tick marks
         'XTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...  %   and tick labels
         'YTick', 1:4, ...
         'YTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...
         'TickLength', [0 0]);
colorbar('southoutside');
colormap(mycolormap);

