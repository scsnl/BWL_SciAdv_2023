clc;
clear all;
close all;


%%% Description %%%%
% This script plot peri-event time-frequency analysis result of neuronal
% GCaMP signals from AI, Cg, PrL, and RSC.

Current_dir=pwd;

%add relevant paths
idcs = strfind(Current_dir,'/');
utils_dir = Current_dir(1:idcs(end-2)-1);
addpath(fullfile(utils_dir,'/','utils'));

data_dir=fullfile(Current_dir, '/', 'spectrogram_analysis');


area='AI';
% area='PrL';
% area='RSC';
% area='Cg';
FigureOutput=0;
time_range=[-2 10];
freq_range=[0.2 1.5]; 
p_threshold=0.02;

range=170:210; %.85-1.05 Hz: Range is all equivalent to all the regions

if exist('X','var')==0
if strcmp(area,'PrL')
    load(fullfile(data_dir,'/','PrL.mat'),'oddball*','trial_ID')
elseif strcmp(area,'Cg')
    load(fullfile(data_dir,'/','Cg.mat'),'oddball*','trial_ID')
elseif strcmp(area,'RSC')
    load(fullfile(data_dir,'/','RSC.mat'),'oddball*','trial_ID')
elseif strcmp(area,'AI')
    load(fullfile(data_dir,'/','AI.mat'),'oddball*','trial_ID')
end
end

if strcmp(area,'AI')
    amplitude_range=[-10 20];
else
    amplitude_range=[-15 15];
end

%%
%%%% universal code to make figure for the publicaion %%%%
% close all;
clear img
% define the figure size and output file name
fig_size=[0 0 12 7]; % bottom left corner of the fig and figure size = [x y width height]
unit='centimeters'; % unit of the figgure & sub-figure size
figure_name='d/dt GCaMP'; % output file name
% pre-set, should be always consistent in one project
alw = 0.75; % AxesLineWidth
fsz = 12;  % Fontsize
lw = 1.5; % LineWidth
msz = 8; % MarkerSize
% creat a new figure
% f=figure('Name',figure_name,'color','white','units',unit,'position',fig_size,'visible','on');
% f.Renderer='Painters';

if exist('X','var')==0
X=[];
if strcmp(area,'AI')
    for id=[1,2,5:length(trial_ID)]
        eval(['X=cat(3,X,',char(trial_ID(id)),'.X(:,:,1:end-3));'])
    end 
else
    for id=1:length(trial_ID)
        eval(['X=cat(3,X,',char(trial_ID(id)),'.X);'])
    end
end
meanX=mean(X,3);
end

y=mean(meanX(range,:)); %original
mean_y=mean(y);
y=y./mean_y;
img(:,:)=mean(X(range,:,:));
clear stderr
stderr=(std(img')./sqrt(size(X,3)-1))./mean_y;
err_p=y+stderr;
err_n=fliplr(y-stderr);
t=([1:401]-200)./10;
t=[t,fliplr(t)];


%%%% universal code to make figure for the publicaion %%%%
%close all; clc;
% define the figure size and output file name
fig_size=[0 13 12 3]; % bottom left corner of the fig and figure size = [x y width height]
unit='centimeters'; % unit of the figgure & sub-figure size
figure_name='Map'; % output file name
% pre-set, should be always consistent in one project
alw = 0.75; % AxesLineWidth
fsz = 12;  % Fontsize
lw = 1.5; % LineWidth
msz = 8; % MarkerSize
% creat a new figure
f3=figure('Name',figure_name,'color','white','units',unit,'position',fig_size,'visible','on');
f3.Renderer='Painters';

if exist('Gmap_seg_percentage','var')==0
%    clear Gmap_seg_percentage
X_mean=mean(X(20:400,:,:),3);
    for i=1:size(X,3)
    temp=X(20:400,:,i);
    Gmap_seg_percentage(:,:,i)=(temp./mean(temp,2))-1;
    %Gmap_seg_percentage(:,:,i)=(temp-mean(X_mean,2))./std(X_mean,[],2);
    end
    
%    clear Gmap_sig_percentage
    for i=1:381
    for j=1:401
    [Gmap_sig_percentage(i,j,1),Gmap_sig_percentage(i,j,2)]=ttest(Gmap_seg_percentage(i,j,:));
    end
    clc
    disp(i)
    end
end
%%

Gmap_seg_percentage_mean=mean(Gmap_seg_percentage,3).*100;
Gmap_nonsig=(Gmap_sig_percentage(:,:,2)>=p_threshold).*Gmap_seg_percentage_mean;
Gmap_nonsig=Gmap_nonsig./max(abs(Gmap_nonsig(:)));
Gmap_nonsig=Gmap_nonsig.*0.5;


if sum(sum(Gmap_sig_percentage(:,:,2)<p_threshold))>0
Gmap_sig_p=(Gmap_sig_percentage(:,:,2)<p_threshold).*(Gmap_seg_percentage_mean>0).*Gmap_seg_percentage_mean;
Gmap_sig_p=(Gmap_sig_p-min(Gmap_sig_p(Gmap_sig_p>0))).*(Gmap_sig_p>0);
Gmap_sig_p=Gmap_sig_p./max(Gmap_sig_p(:))+1;
Gmap_sig_p(Gmap_sig_p==1)=0;
Gmap_sig_n=(Gmap_sig_percentage(:,:,2)<p_threshold).*(Gmap_seg_percentage_mean<0).*Gmap_seg_percentage_mean;
Gmap_sig_n=(Gmap_sig_n-max(Gmap_sig_n(Gmap_sig_n<0))).*(Gmap_sig_n<0);
Gmap_sig_n=-Gmap_sig_n./min(Gmap_sig_n(:))-1;
Gmap_sig_n(Gmap_sig_n==-1)=0;
Gmap_figure=Gmap_sig_p+Gmap_sig_n+Gmap_nonsig;
else
Gmap_figure=Gmap_nonsig;
end
faxis=[0.005:0.005:2];
imagesc([-20.1:0.1:19.9],faxis(20:400),Gmap_figure)
shading flat;
set(gca,'FontSize',12)
xlim(time_range)
set(gca,'xtick',[])
set(gca,'YDir','normal');
ylim(freq_range)
yticks([0.2:0.4:freq_range(2)])
box off
mycolormap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'});
colorbar('southoutside');
caxis([-1.8 1.8])
colormap(mycolormap);
axis off;


if FigureOutput==1
saveas(gcf,['oddball_TFmap_',area],'pdf')
end