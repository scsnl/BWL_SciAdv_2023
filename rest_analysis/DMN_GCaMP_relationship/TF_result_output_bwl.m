function [voting, Gmap_figure, Gmap_sig_p, Gmap_sig_n]=TF_result_output_bwl(TF_subject_means,TF_sig_percentage,p_threshold,cluster_size,pre,post,faxis)

%p_threshold=0.001;
Gmap_seg_percentage_mean=mean(TF_subject_means,3).*100;

sig_mask=(TF_sig_percentage(:,:,2)<p_threshold);
sig_mask = bwareaopen(sig_mask, cluster_size);

Gmap_nonsig=(sig_mask==0).*Gmap_seg_percentage_mean;
Gmap_nonsig=Gmap_nonsig./max(abs(Gmap_nonsig(:)));

Gmap_nonsig=Gmap_nonsig.*0.5; 

if sum(sig_mask(:))>0

    Gmap_sig_p=(sig_mask).*(Gmap_seg_percentage_mean>0).*Gmap_seg_percentage_mean;
    if isempty(Gmap_sig_p(Gmap_sig_p>0))==0
    Gmap_sig_p=(Gmap_sig_p-min(Gmap_sig_p(Gmap_sig_p>0))).*(Gmap_sig_p>0);
    Gmap_sig_p=Gmap_sig_p./max(Gmap_sig_p(:))+1;
    Gmap_sig_p(Gmap_sig_p==1)=0;
    end

    Gmap_sig_n=(sig_mask).*(Gmap_seg_percentage_mean<0).*Gmap_seg_percentage_mean;
    if isempty(Gmap_sig_n(Gmap_sig_n<0))==0
    Gmap_sig_n=(Gmap_sig_n-max(Gmap_sig_n(Gmap_sig_n<0))).*(Gmap_sig_n<0);
    Gmap_sig_n=-Gmap_sig_n./min(Gmap_sig_n(:))-1;
    Gmap_sig_n(Gmap_sig_n==-1)=0;
    end

    Gmap_figure=Gmap_sig_p+Gmap_sig_n+Gmap_nonsig;

else

    Gmap_figure=Gmap_nonsig;
    
end

pcolor((-pre:post)./10,faxis,Gmap_figure)
shading flat;
set(gca,'FontSize',12)
%xlim([-20 20])
set(gca,'xtick',[])
ylim([0.1 1.5])
yticks([0.5:0.5:1.5])
box off

% map_p=[ones(10,1),linspace(0,1,10)',zeros(10,1)];
% map_p_sub=[linspace(0,0.5,10)',zeros(10,1),zeros(10,1)];
% map_n_sub=[zeros(10,1),zeros(10,1),linspace(0.5,0,10)'];
% map_n=[zeros(10,1),linspace(1,0,10)',ones(10,1)];
% map=[map_n;map_n_sub;map_p_sub;map_p];
% colormap(map)
% axis off;
% caxis([-2 2])
% voting=sum(Gmap_sig_p(:,1:100)>0,2);

mycolormap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'});
colorbar('southoutside');
colormap(mycolormap);
axis off;
caxis([-2 2])

% voting=sum(Gmap_sig_p(:,1:100)>0,2);
voting=[-1*sum(Gmap_sig_n(:,1:100)<0,2)+sum(Gmap_sig_p(:,1:100)>0,2)];