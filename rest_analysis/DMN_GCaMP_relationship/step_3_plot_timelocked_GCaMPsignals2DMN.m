close all;


pre=110; post=60; pval=0.05; cluster_size=400;
freq_range=[0.1 1.5];
faxis=linspace(0.01, 1.9, cluster_size);
freq_range_idx=[28:252]; %0.135-1.5Hz

%Correlation between DMN-time locked GCaMP signals
AI_trial_valley=[]; RSC_trial_valley=[]; ACC_trial_valley=[]; mPFC_trial_valley=[]; r_AI_ACC=[]; r_AI_RSC=[]; r_AI_mPFC=[];
AI_trial_peak=[]; RSC_trial_peak=[]; ACC_trial_peak=[]; mPFC_trial_peak=[]; r_AI_ACC=[]; r_AI_RSC=[]; r_AI_mPFC=[];
for run=1:16
    AI_trial_valley(:,run)=mean(squeeze(TF_AI_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);
    RSC_trial_valley(:,run)=mean(squeeze(TF_RSC_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);
    ACC_trial_valley(:,run)=mean(squeeze(TF_ACC_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);
    mPFC_trial_valley(:,run)=mean(squeeze(TF_mPFC_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);

    AI_trial_peak(:,run)=mean(squeeze(TF_AI_GCaMP_peak_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);
    RSC_trial_peak(:,run)=mean(squeeze(TF_RSC_GCaMP_peak_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);
    ACC_trial_peak(:,run)=mean(squeeze(TF_ACC_GCaMP_peak_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);
    mPFC_trial_peak(:,run)=mean(squeeze(TF_mPFC_GCaMP_peak_subject_means_minus_baseline(freq_range_idx,(-pre:post)+615,run)),1);
    
    AI_trial_2valley(:,run)=mean(squeeze(TF_AI_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:0)+615,run)),1);
    RSC_trial_2valley(:,run)=mean(squeeze(TF_RSC_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:0)+615,run)),1);
    ACC_trial_2valley(:,run)=mean(squeeze(TF_ACC_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:0)+615,run)),1);
    mPFC_trial_2valley(:,run)=mean(squeeze(TF_mPFC_GCaMP_valley_subject_means_minus_baseline(freq_range_idx,(-pre:0)+615,run)),1);

    ts=[];
%     ts=[AI_trial_valley(:,run), ACC_trial_valley(:,run), mPFC_trial_valley(:,run), RSC_trial_valley(:,run)];
    ts=[AI_trial_2valley(:,run), ACC_trial_2valley(:,run), mPFC_trial_2valley(:,run), RSC_trial_2valley(:,run)];
    [r_ts(:,:,run), p_ts(:,:,run)]= partialcorr((ts));

end

%% DMN activation
figure;
subplot(5,1,1);
hold off
CBV_seg=-CBV_valley_grp((-pre:10:post)./10+61,:);
plot((-pre:10:post)./10,mean(CBV_seg,2),'k')
hold on
t=[-pre:10:post,fliplr(-pre:10:post)]./10;
err_p=mean(CBV_seg,2)+(std(CBV_seg,[],2));
err_n=flipud(mean(CBV_seg,2)-(std(CBV_seg,[],2)));
fill(t,[err_p;err_n]','b','facealpha',0.3,'linestyle','none')
title('DMN activation in fMRI data')
set(gca,'FontSize',12)
xlim([-pre post]./10)
ylim([-3 3])
box off

subplot(5,1,2);
[voting_AI, Gmap_figure_AI,Gmap_figure_AI_p, Gmap_figure_AI_n]=TF_result_output_bwl(TF_AI_GCaMP_valley_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_AI_GCaMP_valley_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('AI')

subplot(5,1,3);
[voting_ACC, Gmap_figure_ACC, Gmap_figure_ACC_p, Gmap_figure_ACC_n]=TF_result_output_bwl(TF_ACC_GCaMP_valley_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_ACC_GCaMP_valley_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('Cg')

subplot(5,1,4);
[voting_mPFC, Gmap_figure_mPFC, Gmap_figure_mPFC_p, Gmap_figure_mPFC_n]=TF_result_output_bwl(TF_mPFC_GCaMP_valley_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_mPFC_GCaMP_valley_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('PrL')

subplot(5,1,5);
[voting_RSC, Gmap_figure_RSC, Gmap_figure_RSC_p, Gmap_figure_RSC_n]=TF_result_output_bwl(TF_RSC_GCaMP_valley_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_RSC_GCaMP_valley_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('RSC')




%% DMN Deactivation
figure;
subplot(5,1,1);
hold off
CBV_seg=-CBV_peak_grp((-pre:10:post)./10+61,:);
plot((-pre:10:post)./10,mean(CBV_seg,2),'k')
hold on
t=[-pre:10:post,fliplr(-pre:10:post)]./10;
err_p=mean(CBV_seg,2)+(std(CBV_seg,[],2));
err_n=flipud(mean(CBV_seg,2)-(std(CBV_seg,[],2)));
fill(t,[err_p;err_n]','b','facealpha',0.3,'linestyle','none')
title('DMN deactivation in fMRI data')
set(gca,'FontSize',12)
xlim([-pre post]./10)
ylim([-3 3])
box off

subplot(5,1,2);
[voting_AI, Gmap_figure_AI]=TF_result_output_bwl(TF_AI_GCaMP_peak_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_AI_GCaMP_peak_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('AI')

subplot(5,1,3);
[voting_ACC, Gmap_figure_ACC]=TF_result_output_bwl(TF_ACC_GCaMP_peak_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_ACC_GCaMP_peak_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('Cg')

subplot(5,1,4);
[voting_mPFC, Gmap_figure_mPFC]=TF_result_output_bwl(TF_mPFC_GCaMP_peak_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_mPFC_GCaMP_peak_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('PrL')

subplot(5,1,5);
[voting_RSC, Gmap_figure_RSC]=TF_result_output_bwl(TF_RSC_GCaMP_peak_subject_means_minus_baseline(:,(-pre:post)+615,:),TF_RSC_GCaMP_peak_sig_minus_baseline_percentage(:,(-pre:post)+615,:),pval,cluster_size,pre,post,faxis);
ylim(freq_range)
title('RSC')





