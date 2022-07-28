

trial_ID=who('oddball*');

for id=1:length(trial_ID)

disp(id)
    
eval(['t_488=',char(trial_ID(id)),'.t_488;'])
eval(['T=',char(trial_ID(id)),'.T;'])

t_488_zscore=ztrans(t_488(101:end-100)',1:size(t_488(101:end-100),2))';
t_488_high=filter_2sIIR(t_488_zscore,0.005,10,5,'high'); %0.005
eval([char(trial_ID(id)),'.t_488_high=t_488_high;'])

T=T-100;

TF_map = tfa_morlet(t_488_high, 10, 0.005, 2, 0.005);

clear X
pre=200;
post=200;
j=1;
for i=1:length(T)
    if (T(i)-pre) > 0 && (T(i)+post)<size(TF_map,2)
    X(:,:,j)=TF_map(:,(T(i)-pre):(T(i)+post));
    j=j+1;
    end
end

eval([char(trial_ID(id)),'.X=X;'])
end

clearvars -except oddball* trial*