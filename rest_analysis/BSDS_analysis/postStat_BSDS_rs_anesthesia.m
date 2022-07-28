clc
clear all
close all

%%% Description %%%%
% This script load BSDS results and 
% examine spatiotemporal dynamics of latent brain states

Current_dir=pwd;
Model_dir=fullfile(Current_dir,'/', 'model');

%add relevant paths
idcs = strfind(Current_dir,'/');
utils_dir = Current_dir(1:idcs(end-1)-1);
addpath(fullfile(utils_dir,'/','utils'));

% Load model
load(fullfile(Model_dir,'/', 'model_rat_GCaMP_anesthesia_4ROIs_optimal.mat'));
% Load covariance matrices
load(fullfile(Model_dir,'/', 'covMtx_rat_GCaMP_anesthesia_4ROIs_optimal.mat'));


%% Data parameters
TR=0.1; %10Hz
NumbROIs=4;
ROI_names={'AI','Cg', 'PrL', 'RSC'};
num.ROI=NumbROIs;
num.Subj=length(model.temporal_evolution_of_states); % number of subjects.
num.Vol=length(model.temporal_evolution_of_states{1,1}); % length of timeseries
num.State=length(unique(cell2mat(model.temporal_evolution_of_states)));
num.Run=1;


%% Reshape group stats (i.e., temporal evolution of state, posterior probability)
% reshape temporal evolution of state into [subject x run] x [time] & identify dominant states
grpTempEvol=reshape(cell2mat(model.temporal_evolution_of_states),[num.Vol,num.Subj*num.Run])';
dominantStates = unique(grpTempEvol(:));
num.State=length(unique(cell2mat(model.temporal_evolution_of_states)));

% relabel states for ease of computation
grpTempEvol_relabel = zeros(size(grpTempEvol));
for relabel = 1:num.State
    grpTempEvol_relabel(grpTempEvol == dominantStates(relabel)) = relabel;
end

% reshape subjects' posterior probability(pp)
subjStatePP = zeros(num.Subj*num.Run, num.Vol, num.State); % each subject's posterior prob.
idx = 1;
for subj =1:num.Subj
    for run = 1:num.Run
        subjStatePP(idx,:,:) = model.posterior_probabilities{1,subj}((num.Vol*run-(num.Vol-1):num.Vol*run), dominantStates);
        idx = idx+1;
    end
end
grpStatePP=squeeze(mean(subjStatePP,1));% compute group posterior probability by averaging subjects posterior prob.


%% Occupancy rates
OR=[];
for subj=1:num.Subj
    [OR(subj,:),MLT(subj,:),~]=summary_stats_fast(grpTempEvol_relabel(subj,:), 1:num.State);
end


%% Compute State switching matrix
subjStatSwitchProb = zeros(num.Subj, num.State, num.State); % each subject's state switching matrix
for subj=1:num.Subj
    % y-axis:state(t), x-axis:state(t+1)
    SwitchingMatrix = full(sparse(grpTempEvol_relabel(subj, 1:end-1),grpTempEvol_relabel(subj, 2:end),1,num.State,num.State));
    sum_SwitchingMatrix = sum(SwitchingMatrix,2);
    % Normalize the transition matrix
    for state=1:num.State
        SwitchingMatrix(state, :) = SwitchingMatrix(state,:)/sum_SwitchingMatrix(state);
    end
    subjStatSwitchProb(subj,:,:) = SwitchingMatrix;
end
grpStatSwitchProb = squeeze(nanmean(subjStatSwitchProb(1:num.Subj,:,:)));



%% Activation level
Activation=cell(1, 1);
Activation_z_pre=[];
for j= 1:length(dominantStates)
    k = dominantStates(j);
    for subj=1:num.Subj
        Activation{1,j}(subj,:)=model_subjwise.estimated_mean{1,subj}{1,k}';
        Activation_z_pre=[Activation_z_pre; model_subjwise.estimated_mean{1,subj}{1,k}'];
    end
end

%z-score for each region
Activation_z_pre=zscore(Activation_z_pre);
Activation_z=cell(1,1);

for i=1:length(dominantStates)
    Activation_z{1,i}=Activation_z_pre(16*i-15:16*i, :);
end

%% Connectivity analysis: Subject-level
cov_subj=zeros(num.ROI, num.ROI, num.State, num.Subj);
partial_corr_subj=zeros(num.ROI, num.ROI, num.State, num.Subj);
partial_corr_subj_z=zeros(num.ROI, num.ROI, num.State, num.Subj);
partial_corr_mean=zeros(num.ROI, num.ROI, num.State);
partial_corr_sig=zeros(num.ROI, num.ROI, num.State);

for j= 1:length(dominantStates)
    grp_avg_fc=zeros(num.ROI, num.ROI);
    k = dominantStates(j);
    for subj=1:num.Subj
        est_cov = model_subjwise.estimated_covarinace{1,subj}{1,k};
        cov_subj(:,:,j,subj) = est_cov;
        %Partial Correlation
        inv_est_cov = inv(est_cov);
        invD = inv(diag(sqrt(diag(inv_est_cov))));
        partial_corr_subj(:,:,j,subj) = -invD*inv_est_cov*invD;
        partial_corr_subj_z(:,:,j,subj)=1/2.*log((1+partial_corr_subj(:,:,j,subj))./(1-partial_corr_subj(:,:,j,subj)));
    end
    partial_corr_mean(:,:,j)=mean(partial_corr_subj(:,:,j,:),4);
end


p_th=0.01;
p=[];
for k=1:num.State
    for i=1:num.ROI
        for j=1:num.ROI
            [h(i,j,k), p(i,j,k)]=ttest(partial_corr_subj(i,j,k,:));
        end
    end
    p_FDR=FDR(p(:,:,k), p_th);
    if(~isempty(p_FDR))
        partial_corr_sig(:,:,k)=squeeze(partial_corr_mean(:,:,k)).*double(p(:,:,k)<=p_FDR);
        partial_corr_sig(:,:,k)=squeeze(partial_corr_sig(:,:,k)) - diag(diag(squeeze(partial_corr_sig(:,:,k))));
    else
        partial_corr_sig(:,:,k)=zeros(num.ROI, num.ROI);
    end
end



%% Plotting
% reorder the presentation of state based on their occupancy rate.
reorder_seq=[1, 3, 5, 2, 4]; 

%specify colorstrings and colormap for presenting each state
Colorstrings = {'#281B24','#3E3842','#74858F','#A4BDBE','#BDD6DA'};
Colormap=[40 27 36;...
          62 56 66;...
          116 133 143;...
          164 189 190;...
          189 214 218]./255;


%% Plot time-varying posterior probability & temporal evolution
figure;
subplot(2,1,1);
plot(grpStatePP);
for i = 1:num.State
    hold on;
    plot(grpStatePP(:,i), 'Color', Colorstrings{i}, 'LineWidth', 2);
end
ylim([0,1]); xlim([0, num.Vol]); xlabel('Time(s)', 'FontSize', 12); ylabel('posterior probability', 'FontSize', 12);
title('Posterior probability of latent brain states');

subplot(2,1,2);
hAxes = gca;
imagesc(hAxes, grpTempEvol_relabel);
colormap( hAxes, Colormap);
title('Temporal Evolution of latent brain states');


%% Plot occupancy rates
%change order
OR_reorder=OR(:, reorder_seq);
figure;
b=bar(1:num.State, mean(OR_reorder)*100, 'FaceColor', 'flat');
for k=1:num.State
    b.CData(k,:) = Colormap(k,:);
end
hold on; errorbar(1:num.State, mean(OR_reorder).*100, (std(OR_reorder)/sqrt(length(OR_reorder))).*100, 'LineStyle', 'none', 'Color', 'k');  % plot standard error.
set(gca, 'xticklabel', {'S1', 'S2', 'S3', 'S4', 'S5'});
ylim([0 40]);
ylabel('Percentage');
title('Occupancy rate in anesthetized rats');

h_OR=[]; p_OR=[];
for i=1:num.State
    for j=1:num.State
        [h_OR(i,j), p_OR(i,j)]=ttest(OR_reorder(:,i),OR_reorder(:,j));
    end
end
p_FDR=0.01;
p_FDR=FDR(p_OR, p_th);
p_OR_FDR_corrected=double(p_OR<p_FDR);


%% Plot state switching matrix
%change order
grpStatSwitchProb_reorder=grpStatSwitchProb(reorder_seq, reorder_seq);

figure;
imagesc(grpStatSwitchProb_reorder);
ylabel('state(t)');
xlabel('state(t+1)');
title('State transition probability');
caxis([0 1]);
textStrings = num2str(grpStatSwitchProb_reorder(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:num.State);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
% Choose white or black for the text color of the strings so they can be easily seen over the background color
textColors = repmat(grpStatSwitchProb_reorder(:) > midValue, 1, 3);  
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
set(gca, 'XTick', 1:num.State, ...                             % Change the axes tick marks
         'XTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...  %   and tick labels
         'YTick', 1:num.State, ...
         'YTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...
         'TickLength', [0 0]);
mycolormap = customcolormap([0 0.5 1], {'#9d0142','#f66e45','#ffffbb'});
colorbar('southoutside');
colormap(mycolormap);


%% Plot activation
Activation_z_reorder=Activation_z(reorder_seq);

mean_Activation_z_reorder=[];
for i=1:num.State
    mean_Activation_z_reorder(:,i)=mean(Activation_z_reorder{1,i});
end

figure;
for i=1:num.ROI
    hold on;
    plot(mean_Activation_z_reorder(i,:));
end
legend({'AI', 'Cg', 'PrL', 'RSC'});
title('Normalized GCaMP intensity across states')
ylim([-1.5 1.5]);

%Reorder once more to show cyclic state transition
mean_Activation_z_reorder2=mean_Activation_z_reorder(:, [1,5,3,2,4]);
mean_Activation_z_reorder2_plot=repmat(mean_Activation_z_reorder2,1,3);

figure;
for i=1:num.ROI
    hold on;
    plot(mean_Activation_z_reorder2_plot(i,:), '-o');
end
legend({'AI', 'Cg', 'PrL', 'RSC'});
title('Cyclic chcanges of GCaMP intensity through stat evolution')
ylim([-1.5 1.5]);


%% plot partial correlation of each state
partial_corr_subj_reorder_z=partial_corr_subj_z(:,:,reorder_seq,:);

% Plot connectivity changes according to cyclic state transition
partial_corr_subj_reorder2_z=partial_corr_subj_reorder_z(:,:,[1,5,3,2,4],:);
partial_corr_subj_reorder2_z=squeeze(mean(partial_corr_subj_reorder2_z,4));
partial_corr_subj_reorder2_z=repmat(partial_corr_subj_reorder2_z,1,1,3);

figure;
for state=1:5
subplot(1,5,state);
partialcorr_subj=squeeze([partial_corr_subj_reorder_z(1,2,state,:), partial_corr_subj_reorder_z(1,3,state,:), partial_corr_subj_reorder_z(1,4,state,:),...
                          partial_corr_subj_reorder_z(2,3,state,:), partial_corr_subj_reorder_z(2,4,state,:), partial_corr_subj_reorder_z(3,4,state,:)])'; 
groupIdx=[];
for group=1:6
    groupIdx=[groupIdx; group.*ones(num.Subj,1)];
end
plot_box_scatter(partialcorr_subj, groupIdx)
set(gca, 'xticklabel', {'AI-Cg', 'AI-PrL', 'AI-RSC', 'Cg-PrL', 'Cg-RSC', 'PrL-RSC'});
title(sprintf('Connectivity of State %d', state));
ylim([-0.4 1]);
end
sgtitle('Functional connectivity pattern of latent brain states')

figure;
hold on;
plot(squeeze(partial_corr_subj_reorder2_z(1,2,:)), '-o');
hold on;
plot(squeeze(partial_corr_subj_reorder2_z(1,3,:)), '-o');
hold on;
plot(squeeze(partial_corr_subj_reorder2_z(1,4,:)), '-o');
hold on;
plot(squeeze(partial_corr_subj_reorder2_z(2,3,:)), '-o');
hold on;
plot(squeeze(partial_corr_subj_reorder2_z(2,4,:)), '-o');
hold on;
plot(squeeze(partial_corr_subj_reorder2_z(3,4,:)), '-o');
legend({'AI-Cg', 'AI-PrL', 'AI-RSC', 'Cg-PrL', 'Cg-RSC', 'PrL-RSC'});
title('Cyclic changes of inter-regional functional connectivity');



