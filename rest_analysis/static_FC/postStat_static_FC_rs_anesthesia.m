clc
clear all 
close all

%%% Description %%%%
% This script computes static, time-averaged functional connectivity (partial correlation)
% between regions during resting state in anesthetized rats.

Current_dir=pwd;

%add relevant paths
idcs = strfind(Current_dir,'/');
addpath(fullfile(Current_dir(1:idcs(end-1)-1),'/','utils'));

% Load data. Order: AI, Cg, PrL, RSC
Data_dir=fullfile(Current_dir(1:idcs(end)-1),'/', 'data');
load(fullfile(Data_dir, '/', 'GCaMP_anesthesia_resting.mat'));

%% Data parameters
TR=0.1; %10Hz
NumbROIs=4;
num.ROI=NumbROIs;
num.Subj=length(data);
 
%% Compute partial correlation between regions.
for subj=1:num.Subj
    [r_pc(:,:,subj), ~]=partialcorr(data{1,subj}');
end

for i=1:num.ROI
    for j=1:num.ROI
        [h_pc(i,j), p_pc(i,j)]=ttest(r_pc(i,j,:));
    end
end

r_pc_avg=mean(r_pc,3);

p_th=0.01;
p_FDR_pc=FDR(p_pc, p_th);
if(~isempty(p_FDR_pc))
    r_pc_sig=r_pc_avg.*double(p_pc<=p_FDR_pc);
    r_pc_sig=r_pc_sig - diag(diag(r_pc_sig));
else
    r_pc_sig=zeros(num.ROI, num.ROI);
end

groupIdx=[];
for group=1:6
    groupIdx=[groupIdx; group.*ones(subj,1)];
end
figure;
partialcorr_subj=squeeze([r_pc(1,2,:), r_pc(1,3,:), r_pc(1,4,:), r_pc(2,3,:), r_pc(2,4,:), r_pc(3,4,:)])'; 
plot_box_scatter(partialcorr_subj, groupIdx)
set(gca, 'xticklabel', {'AI-Cg', 'AI-PrL', 'AI-RSC', 'Cg-PrL', 'Cg-RSC', 'PrL-RSC'});
title('Partial correlation between regions');

dist=sqrt(2.*(1-abs(r_pc_avg)));
dist=dist - diag(diag(dist));
z=linkage(dist);
figure;
dendrogram(z);
title('Hierarchical organization of functional connectivity between putative rodent DMN nodes and AI')







