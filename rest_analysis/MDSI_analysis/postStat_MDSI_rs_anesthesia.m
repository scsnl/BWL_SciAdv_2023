clc
clear all
close all

%%% Description %%%%
% This script load MDSI inference results and do statistical analysis to
% identify significant causal interactions during resting state in anesthetized rats.

Current_dir=pwd;
Model_Dir=fullfile(Current_dir,'/', 'model');

state='anesthesia';
NumbSamples=16;
    
NumbROI=4;
for trial=1:NumbSamples
    load(fullfile(Model_Dir, sprintf('Conn_GCaMP_%s_resting_trial_%d.mat', state, trial)));
    AC_trial(:,:,trial)=A;
end
        
h=[]; p=[];
AC_mean=[];
AC_mean(:,:)=mean(AC_trial,3);

for i=1:NumbROI
    for j=1:NumbROI
        [h(i,j), p(i,j)]=ttest(AC_trial(i,j,:));
    end
end

p_th=0.05;

offdiag_idx=[1:NumbROI^2];
offdiag_idx(find(eye(NumbROI, NumbROI)))=[];
p_FDR=FDR(p(offdiag_idx), p_th);
if(~isempty(p_FDR))
    AC_sig=AC_mean.*double(p<=p_FDR);
    AC_sig=AC_sig - diag(diag(AC_sig));
else
    AC_sig=zeros(NumbROI, NumbROI);
end

% Plot causal interaction of each task conditions
figure;
imagesc(AC_sig);
ylabel('To');
xlabel('From');
title(sprintf('Causal interactions in %s (FDR corrected)',state));
caxis([-0.1 0.1]);
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
set(gca, 'XTick', 1:4, ...                             % Change the axes tick marks
         'XTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...  %   and tick labels
         'YTick', 1:4, ...
         'YTickLabel', {'AI', 'Cg', 'PrL', 'RSC'}, ...
         'TickLength', [0 0]);
colorbar('southoutside');
colormap(mycolormap);



