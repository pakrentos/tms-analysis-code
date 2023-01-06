load Results_sources_2nd_Day\4-8Hz\sources_Pre_RH_Im1_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta_RSNorm.mat;
load Results_sources_2nd_Day\4-8Hz\sources_Post_RH_Im1_foi=6Hz_df=2Hz_5-5.5s_eLoreta_RSNorm.mat;
pre1 = sourcepre_Norm;
post1 = sourcepost_Norm;

load Results_sources_sham_Day\4-8Hz\sources_Pre_RH_Im1_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta_RSNorm.mat;
load Results_sources_sham_Day\4-8Hz\sources_Post_RH_Im1_foi=6Hz_df=2Hz_5-5.5s_eLoreta_RSNorm.mat;
pre1(21:35) = sourcepre_Norm;
post1(21:35) = sourcepost_Norm;
 
sourcepre_Norm = pre1;
sourcepost_Norm = post1;

stat.posmask         = (stat.negclusterslabelmat<2)&(stat.negclusterslabelmat>0);

%stat.mask2         = (stat.stat>0.7*max(stat.stat));
% stat.mask2         = (stat.stat<0.7*min(stat.stat));
% stat.mask = stat.mask2.*stat.posmask;

subj_set=[1:9 11:35];

kpow=find(stat.posmask==1);

for subj=subj_set
    pow_avg(subj,1) = mean(sourcepre_Norm(subj).avg.pow(kpow));
    pow_avg(subj,2) = mean(sourcepost_Norm(subj).avg.pow(kpow));
end;

csvwrite('ClustPow_new\pow_in_clust_Sham_and_TMS_Im1_theta_f0=6z_df=2Hz_Pre_Post_CA=0.02.csv', pow_avg);

%%
for n=1:2
    pow_grAvg(n)=0;
    for subj=subj_set
        pow_grAvg(n)=pow_grAvg(n)+pow_avg(subj,n);
    end;
    pow_grAvg(n)=pow_grAvg(n)/29;
end;

%%
figure; plot(pow_grAvg)
