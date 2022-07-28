baseline_range=1:1201; % include all. Chose to go with this baselinecovering entire signal

% mPFC peak
TF_mPFC_GCaMP_peak_grp=TF_mPFC_GCaMP_peak_grp./mean(TF_mPFC_GCaMP_peak_grp,2);
start=1;
for run=1:16
TF_mPFC_GCaMP_peak_subject_means(:,:,run)=mean(TF_mPFC_GCaMP_peak_grp(:,:,start:start+PeakNums(run)-1),3);
if run<32
start=start+PeakNums(run);
end
end

TF_mPFC_GCaMP_peak_subject_means_minus_baseline=TF_mPFC_GCaMP_peak_subject_means./mean(TF_mPFC_GCaMP_peak_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_mPFC_GCaMP_peak_sig_minus_baseline_percentage(i,j,1),TF_mPFC_GCaMP_peak_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_mPFC_GCaMP_peak_subject_means_minus_baseline(i,j,:));
end
clc
disp('mPFC')
disp(i)
end

% ACC peak
TF_ACC_GCaMP_peak_grp=TF_ACC_GCaMP_peak_grp./mean(TF_ACC_GCaMP_peak_grp,2);
start=1;
for run=1:16
TF_ACC_GCaMP_peak_subject_means(:,:,run)=mean(TF_ACC_GCaMP_peak_grp(:,:,start:start+PeakNums(run)-1),3);
if run<32
start=start+PeakNums(run);
end
end

TF_ACC_GCaMP_peak_subject_means_minus_baseline=TF_ACC_GCaMP_peak_subject_means./mean(TF_ACC_GCaMP_peak_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_ACC_GCaMP_peak_sig_minus_baseline_percentage(i,j,1),TF_ACC_GCaMP_peak_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_ACC_GCaMP_peak_subject_means_minus_baseline(i,j,:));
end
clc
disp('ACC')
disp(i)
end

% RSC peak
TF_RSC_GCaMP_peak_grp=TF_RSC_GCaMP_peak_grp./mean(TF_RSC_GCaMP_peak_grp,2);
start=1;
for run=1:16
TF_RSC_GCaMP_peak_subject_means(:,:,run)=mean(TF_RSC_GCaMP_peak_grp(:,:,start:start+PeakNums(run)-1),3);
if run<32
start=start+PeakNums(run);
end
end

TF_RSC_GCaMP_peak_subject_means_minus_baseline=TF_RSC_GCaMP_peak_subject_means./mean(TF_RSC_GCaMP_peak_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_RSC_GCaMP_peak_sig_minus_baseline_percentage(i,j,1),TF_RSC_GCaMP_peak_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_RSC_GCaMP_peak_subject_means_minus_baseline(i,j,:));
end
clc
disp('RSC')
disp(i)
end

% AI peak
TF_AI_GCaMP_peak_grp=TF_AI_GCaMP_peak_grp./mean(TF_AI_GCaMP_peak_grp,2);
start=1;
for run=1:16
TF_AI_GCaMP_peak_subject_means(:,:,run)=mean(TF_AI_GCaMP_peak_grp(:,:,start:start+PeakNums(run)-1),3);
if run<32
start=start+PeakNums(run);
end
end

TF_AI_GCaMP_peak_subject_means_minus_baseline=TF_AI_GCaMP_peak_subject_means./mean(TF_AI_GCaMP_peak_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_AI_GCaMP_peak_sig_minus_baseline_percentage(i,j,1),TF_AI_GCaMP_peak_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_AI_GCaMP_peak_subject_means_minus_baseline(i,j,:));
end
clc
disp('AI')
disp(i)
end
%%
% mPFC valley
TF_mPFC_GCaMP_valley_grp=TF_mPFC_GCaMP_valley_grp./mean(TF_mPFC_GCaMP_valley_grp,2);
start=1;
for run=1:16
TF_mPFC_GCaMP_valley_subject_means(:,:,run)=mean(TF_mPFC_GCaMP_valley_grp(:,:,start:start+ValleyNums(run)-1),3);
if run<32
start=start+ValleyNums(run);
end
end

TF_mPFC_GCaMP_valley_subject_means_minus_baseline=TF_mPFC_GCaMP_valley_subject_means./mean(TF_mPFC_GCaMP_valley_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_mPFC_GCaMP_valley_sig_minus_baseline_percentage(i,j,1),TF_mPFC_GCaMP_valley_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_mPFC_GCaMP_valley_subject_means_minus_baseline(i,j,:));
end
clc
disp('mPFC')
disp(i)
end

% ACC valley
TF_ACC_GCaMP_valley_grp=TF_ACC_GCaMP_valley_grp./mean(TF_ACC_GCaMP_valley_grp,2);
start=1;
for run=1:16
TF_ACC_GCaMP_valley_subject_means(:,:,run)=mean(TF_ACC_GCaMP_valley_grp(:,:,start:start+ValleyNums(run)-1),3);
if run<32
start=start+ValleyNums(run);
end
end

TF_ACC_GCaMP_valley_subject_means_minus_baseline=TF_ACC_GCaMP_valley_subject_means./mean(TF_ACC_GCaMP_valley_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_ACC_GCaMP_valley_sig_minus_basCBV_valley_grp_subject_meanseline_percentage(i,j,1),TF_ACC_GCaMP_valley_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_ACC_GCaMP_valley_subject_means_minus_baseline(i,j,:));
end
clc
disp('ACC')
disp(i)
end

%RSC valley
TF_RSC_GCaMP_valley_grp=TF_RSC_GCaMP_valley_grp./mean(TF_RSC_GCaMP_valley_grp,2);
start=1;
for run=1:16
TF_RSC_GCaMP_valley_subject_means(:,:,run)=mean(TF_RSC_GCaMP_valley_grp(:,:,start:start+ValleyNums(run)-1),3);
if run<32
start=start+ValleyNums(run);
end
end

TF_RSC_GCaMP_valley_subject_means_minus_baseline=TF_RSC_GCaMP_valley_subject_means./mean(TF_RSC_GCaMP_valley_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_RSC_GCaMP_valley_sig_minus_baseline_percentage(i,j,1),TF_RSC_GCaMP_valley_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_RSC_GCaMP_valley_subject_means_minus_baseline(i,j,:));
end
clc
disp('RSC')
disp(i)
end

% AI valley
TF_AI_GCaMP_valley_grp=TF_AI_GCaMP_valley_grp./mean(TF_AI_GCaMP_valley_grp,2);
start=1;
for run=1:16
TF_AI_GCaMP_valley_subject_means(:,:,run)=mean(TF_AI_GCaMP_valley_grp(:,:,start:start+ValleyNums(run)-1),3);
if run<32
start=start+ValleyNums(run);
end
end

TF_AI_GCaMP_valley_subject_means_minus_baseline=TF_AI_GCaMP_valley_subject_means./mean(TF_AI_GCaMP_valley_subject_means(:,baseline_range,:),2)-1;
for i=1:400
for j=1:1201
[TF_AI_GCaMP_valley_sig_minus_baseline_percentage(i,j,1),TF_AI_GCaMP_valley_sig_minus_baseline_percentage(i,j,2)]=ttest(TF_AI_GCaMP_valley_subject_means_minus_baseline(i,j,:));
end
clc
disp('AI')
disp(i)
end