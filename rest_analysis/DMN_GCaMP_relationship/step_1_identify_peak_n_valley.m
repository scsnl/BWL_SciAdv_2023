clc; clear all; close all;

Current_dir=pwd;


%add relevant paths
idcs = strfind(Current_dir,'/');
addpath(fullfile(Current_dir(1:idcs(end-1)-1),'/','utils'));

% Load data. Order: AI, Cg, PrL, RSC
Data_dir=fullfile(Current_dir(1:idcs(end)-1),'/', 'data','/', 'spectrogram_analysis');

CBV_peak_shift=0;
CBV_valley_shift=0;
N=600; %total fMRI frame
pre_for_epi=60; post_for_epi=60; % seg range
% D=dir('DMN*10.1D'); % DMN*08.1D or DMN*10.1D
D=dir(fullfile(Data_dir, 'DMN*10.1D'));

TF_mPFC_GCaMP_peak_grp=[];
TF_ACC_GCaMP_peak_grp=[];
TF_RSC_GCaMP_peak_grp=[];
TF_AI_GCaMP_peak_grp=[];
mPFC_GCaMP_peak_grp=[];
ACC_GCaMP_peak_grp=[];
RSC_GCaMP_peak_grp=[];
AI_GCaMP_peak_grp=[];
mPFC_Rho_peak_grp=[];
ACC_Rho_peak_grp=[];
RSC_Rho_peak_grp=[];
AI_Rho_peak_grp=[];
CBV_peak_grp=[];
CBV_peak_index_grp=[];

TF_mPFC_GCaMP_valley_grp=[];
TF_ACC_GCaMP_valley_grp=[];
TF_RSC_GCaMP_valley_grp=[];
TF_AI_GCaMP_valley_grp=[];
mPFC_GCaMP_valley_grp=[];
ACC_GCaMP_valley_grp=[];
RSC_GCaMP_valley_grp=[];
AI_GCaMP_valley_grp=[];
mPFC_Rho_valley_grp=[];
ACC_Rho_valley_grp=[];
RSC_Rho_valley_grp=[];
AI_Rho_valley_grp=[];
CBV_valley_grp=[];
CBV_valley_index_grp=[];

PeakNums=zeros(1,length(D));
ValleyNums=zeros(1,length(D));

bwl=load(fullfile(Data_dir, '/','DMN_TimeSeries&PeaksEvents.mat'), 'CBV_grp', 'CBV_peak_index_grp', 'CBV_valley_index_grp');


for run=1:length(D)
disp([num2str(length(D)),'-',num2str(run)])
load(fullfile(Data_dir, '/', 'ALL_EPI_GCaMP_and_Rho.mat'),['*',D(run).name(1:9),'0',D(run).name(10),'*'])

%%
% DMN_amplitude=textread(D(run).name);
DMN_amplitude=textread(fullfile(Data_dir, '/', D(run).name));

eval(['TF_mPFC_GCaMP_full=',char(who('TF_DMN*mPFC*GCaMP')),';'])
eval(['TF_ACC_GCaMP_full=',char(who('TF_DMN*ACC*GCaMP')),';'])
eval(['TF_RSC_GCaMP_full=',char(who('TF_DMN*RSC*GCaMP')),';'])
eval(['TF_AI_GCaMP_full=',char(who('TF_DMN*AI*GCaMP')),';'])

eval(['mPFC_GCaMP_full=',char(who('DMN*mPFC*_GCaMP')),';'])
eval(['ACC_GCaMP_full=',char(who('DMN*ACC*_GCaMP')),';'])
eval(['RSC_GCaMP_full=',char(who('DMN*RSC*_GCaMP')),';'])
eval(['AI_GCaMP_full=',char(who('DMN*AI*_GCaMP')),';'])

eval(['mPFC_Rho_full=',char(who('DMN*mPFC*_Rho')),';'])
eval(['ACC_Rho_full=',char(who('DMN*ACC*_Rho')),';'])
eval(['RSC_Rho_full=',char(who('DMN*RSC*_Rho')),';'])
eval(['AI_Rho_full=',char(who('DMN*AI*_Rho')),';'])

%% Plot full GCaMP signal
% freq_range_idx=[28:252]; %0.135-1.2Hz
% figure;
% subplot(1,4,1);
% plot(zscore(mean(TF_AI_GCaMP_full(freq_range_idx,:)))); title('AI');
% subplot(1,4,2);
% plot(zscore(mean(TF_ACC_GCaMP_full(freq_range_idx,:)))); title('Cg');
% subplot(1,4,3);
% plot(zscore(mean(TF_mPFC_GCaMP_full(freq_range_idx,:)))); title('PrL');
% subplot(1,4,4);
% plot(zscore(mean(TF_RSC_GCaMP_full(freq_range_idx,:)))); title('RSC');


CBV=bwl.CBV_grp{run};
CBV_peak_index=bwl.CBV_peak_index_grp{run};
CBV_valley_index=bwl.CBV_valley_index_grp{run};

% % CBV=filter_2sIIR(DMN_amplitude',[0.01 0.1],1,3,'bandpass')';
% % dCBV=CBV(2:end)-CBV(1:end-1);
% % ddCBV=dCBV(2:end)-dCBV(1:end-1);
% % CBV=(CBV-mean(CBV))./std(CBV);

% %%%%% HIGH
% CBV_peak=(abs(dCBV(1:end-1))<0.05).*(ddCBV<0).*(CBV(1:end-2)>1.63);
% CBV_peak_left=(CBV_peak(1:end-1)-CBV_peak(2:end))<0;
% CBV_peak_left_index=find(CBV_peak_left==1);
% CBV_peak_right=(CBV_peak(1:end-1)-CBV_peak(2:end))>0;
% CBV_peak_right_index=find(CBV_peak_right==1);
% CBV_peak_index=zeros(1,length(CBV_peak_left_index));
% for peak_i=1:min([length(CBV_peak_left_index),length(CBV_peak_right_index)])
% CBV_peak_index(peak_i)=find(CBV(CBV_peak_left_index(peak_i):CBV_peak_right_index(peak_i))==max(CBV(CBV_peak_left_index(peak_i):CBV_peak_right_index(peak_i))))+CBV_peak_left_index(peak_i)-1;
% end
% CBV_peak_index=CBV_peak_index+CBV_peak_shift; % CBV_peak_shift
% CBV_peak_index(CBV_peak_index<=pre_for_epi)=[];
% CBV_peak_index(CBV_peak_index>=(N-post_for_epi))=[];
% CBV_peak_index=CBV_peak_index.*10;
% 
% %%%%% LOW
% CBV_valley=(abs(dCBV(1:end-1))<0.05).*(ddCBV>0).*(CBV(1:end-2)<-1.63);
% CBV_valley_left=(CBV_valley(1:end-1)-CBV_valley(2:end))<0;
% CBV_valley_left_index=find(CBV_valley_left==1);
% CBV_valley_right=(CBV_valley(1:end-1)-CBV_valley(2:end))>0;
% CBV_valley_right_index=find(CBV_valley_right==1);
% CBV_valley_index=zeros(1,length(CBV_valley_left_index));
% for peak_i=1:min([length(CBV_valley_left_index),length(CBV_valley_right_index)])
%     CBV_valley_index(peak_i)=find(CBV(CBV_valley_left_index(peak_i):CBV_valley_right_index(peak_i))==min(CBV(CBV_valley_left_index(peak_i):CBV_valley_right_index(peak_i))))+CBV_valley_left_index(peak_i)-1;
% end
% CBV_valley_index=CBV_valley_index+CBV_valley_shift; % CBV_valley_shift
% CBV_valley_index(CBV_valley_index<=pre_for_epi)=[];
% CBV_valley_index(CBV_valley_index>=(N-post_for_epi))=[];
% CBV_valley_index=CBV_valley_index.*10;


%
pre=pre_for_epi.*10;
post=post_for_epi.*10;

for peak_i=1:length(CBV_peak_index)
TF_mPFC_GCaMP_peak_individual(:,:,peak_i)=TF_mPFC_GCaMP_full(:,CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
TF_ACC_GCaMP_peak_individual(:,:,peak_i)=TF_ACC_GCaMP_full(:,CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
TF_RSC_GCaMP_peak_individual(:,:,peak_i)=TF_RSC_GCaMP_full(:,CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
TF_AI_GCaMP_peak_individual(:,:,peak_i)=TF_AI_GCaMP_full(:,CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
mPFC_GCaMP_peak_individual(:,peak_i)=mPFC_GCaMP_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
ACC_GCaMP_peak_individual(:,peak_i)=ACC_GCaMP_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
RSC_GCaMP_peak_individual(:,peak_i)=RSC_GCaMP_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
AI_GCaMP_peak_individual(:,peak_i)=AI_GCaMP_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
mPFC_Rho_peak_individual(:,peak_i)=mPFC_Rho_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
ACC_Rho_peak_individual(:,peak_i)=ACC_Rho_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
RSC_Rho_peak_individual(:,peak_i)=RSC_Rho_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
AI_Rho_peak_individual(:,peak_i)=AI_Rho_full(CBV_peak_index(peak_i)-pre:CBV_peak_index(peak_i)+post);
CBV_peak_individual(:,peak_i)=CBV((CBV_peak_index(peak_i)-pre)./10:(CBV_peak_index(peak_i)+post)./10);
end

for valley_i=1:length(CBV_valley_index)
TF_mPFC_GCaMP_valley_individual(:,:,valley_i)=TF_mPFC_GCaMP_full(:,CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
TF_ACC_GCaMP_valley_individual(:,:,valley_i)=TF_ACC_GCaMP_full(:,CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
TF_RSC_GCaMP_valley_individual(:,:,valley_i)=TF_RSC_GCaMP_full(:,CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
TF_AI_GCaMP_valley_individual(:,:,valley_i)=TF_AI_GCaMP_full(:,CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
mPFC_GCaMP_valley_individual(:,valley_i)=mPFC_GCaMP_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
ACC_GCaMP_valley_individual(:,valley_i)=ACC_GCaMP_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
RSC_GCaMP_valley_individual(:,valley_i)=RSC_GCaMP_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
AI_GCaMP_valley_individual(:,valley_i)=AI_GCaMP_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
mPFC_Rho_valley_individual(:,valley_i)=mPFC_Rho_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
ACC_Rho_valley_individual(:,valley_i)=ACC_Rho_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
RSC_Rho_valley_individual(:,valley_i)=RSC_Rho_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
AI_Rho_valley_individual(:,valley_i)=AI_Rho_full(CBV_valley_index(valley_i)-pre:CBV_valley_index(valley_i)+post);
CBV_valley_individual(:,valley_i)=CBV((CBV_valley_index(valley_i)-pre)./10:(CBV_valley_index(valley_i)+post)./10);
end

%%
TF_mPFC_GCaMP_peak_grp=cat(3,TF_mPFC_GCaMP_peak_grp,TF_mPFC_GCaMP_peak_individual);
TF_ACC_GCaMP_peak_grp=cat(3,TF_ACC_GCaMP_peak_grp,TF_ACC_GCaMP_peak_individual);
TF_RSC_GCaMP_peak_grp=cat(3,TF_RSC_GCaMP_peak_grp,TF_RSC_GCaMP_peak_individual);
TF_AI_GCaMP_peak_grp=cat(3,TF_AI_GCaMP_peak_grp,TF_AI_GCaMP_peak_individual);
mPFC_GCaMP_peak_grp=cat(2,mPFC_GCaMP_peak_grp,mPFC_GCaMP_peak_individual);
ACC_GCaMP_peak_grp=cat(2,ACC_GCaMP_peak_grp,ACC_GCaMP_peak_individual);
RSC_GCaMP_peak_grp=cat(2,RSC_GCaMP_peak_grp,RSC_GCaMP_peak_individual);
AI_GCaMP_peak_grp=cat(2,AI_GCaMP_peak_grp,AI_GCaMP_peak_individual);
mPFC_Rho_peak_grp=cat(2,mPFC_Rho_peak_grp,mPFC_Rho_peak_individual);
ACC_Rho_peak_grp=cat(2,ACC_Rho_peak_grp,ACC_Rho_peak_individual);
RSC_Rho_peak_grp=cat(2,RSC_Rho_peak_grp,RSC_Rho_peak_individual);
AI_Rho_peak_grp=cat(2,AI_Rho_peak_grp,AI_Rho_peak_individual);
CBV_peak_grp=cat(2,CBV_peak_grp,CBV_peak_individual);

TF_mPFC_GCaMP_valley_grp=cat(3,TF_mPFC_GCaMP_valley_grp,TF_mPFC_GCaMP_valley_individual);
TF_ACC_GCaMP_valley_grp=cat(3,TF_ACC_GCaMP_valley_grp,TF_ACC_GCaMP_valley_individual);
TF_RSC_GCaMP_valley_grp=cat(3,TF_RSC_GCaMP_valley_grp,TF_RSC_GCaMP_valley_individual);
TF_AI_GCaMP_valley_grp=cat(3,TF_AI_GCaMP_valley_grp,TF_AI_GCaMP_valley_individual);
mPFC_GCaMP_valley_grp=cat(2,mPFC_GCaMP_valley_grp,mPFC_GCaMP_valley_individual);
ACC_GCaMP_valley_grp=cat(2,ACC_GCaMP_valley_grp,ACC_GCaMP_valley_individual);
RSC_GCaMP_valley_grp=cat(2,RSC_GCaMP_valley_grp,RSC_GCaMP_valley_individual);
AI_GCaMP_valley_grp=cat(2,AI_GCaMP_valley_grp,AI_GCaMP_valley_individual);
mPFC_Rho_valley_grp=cat(2,mPFC_Rho_valley_grp,mPFC_Rho_valley_individual);
ACC_Rho_valley_grp=cat(2,ACC_Rho_valley_grp,ACC_Rho_valley_individual);
RSC_Rho_valley_grp=cat(2,RSC_Rho_valley_grp,RSC_Rho_valley_individual);
AI_Rho_valley_grp=cat(2,AI_Rho_valley_grp,AI_Rho_valley_individual);
CBV_valley_grp=cat(2,CBV_valley_grp,CBV_valley_individual);

%%
PeakNums(run)=length(CBV_peak_index);
ValleyNums(run)=length(CBV_valley_index);

CBV_peak_index_grp{run}=CBV_peak_index;
CBV_valley_index_grp{run}=CBV_valley_index;

mPFC_GCaMP_grp{run}=mPFC_GCaMP_full;
ACC_GCaMP_grp{run}=ACC_GCaMP_full;
RSC_GCaMP_grp{run}=RSC_GCaMP_full;
AI_GCaMP_grp{run}=AI_GCaMP_full;

mPFC_Rho_grp{run}=mPFC_Rho_full;
ACC_Rho_grp{run}=ACC_Rho_full;
RSC_Rho_grp{run}=RSC_Rho_full;
AI_Rho_grp{run}=AI_Rho_full;

CBV_grp{run}=CBV;

%%
clear(['*',D(run).name(1:9),'0',D(run).name(10),'*'])
clear *_individual
end

clearvars -except *grp PeakNums ValleyNums D

%% plot low-passed CBV and peak/valley labels
% for run=1:16
% subplot(4,4,run)
% hold off
% plot(CBV_grp{1, run})
% hold on
% plot(CBV_valley_index_grp{1, run}./10,CBV_grp{1, run}(CBV_valley_index_grp{1, run}./10),'o')
% plot(CBV_peak_index_grp{1, run}./10,CBV_grp{1, run}(CBV_peak_index_grp{1, run}./10),'o')
% plot([0 600],[1.63 1.63],':')
% plot([0 600],[-1.63 -1.63],':')
% ylim([-5 5])
% title([D(run).name(1:5),' ',D(run).name(7:10)])
% end

for run=1:16
subplot(4,4,run)
hold off
plot(-1*CBV_grp{1, run})
hold on
plot(CBV_valley_index_grp{1, run}./10, -1*CBV_grp{1, run}(CBV_valley_index_grp{1, run}./10),'o')
plot(CBV_peak_index_grp{1, run}./10, -1*CBV_grp{1, run}(CBV_peak_index_grp{1, run}./10),'o')
plot([0 600],[1.63 1.63],':')
plot([0 600],[-1.63 -1.63],':')
ylim([-5 5])
title([D(run).name(1:5),' ',D(run).name(7:10), '-Fliped'])
end





