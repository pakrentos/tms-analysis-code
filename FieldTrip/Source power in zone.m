global ft_default
ft_default.trackusage = 'no';
ft_defaults;

atlas = ft_read_atlas('ROI_MNI_V4.nii');
atlas_t= ft_convert_coordsys(atlas, 'ctf');
atlas_t = ft_convert_units(atlas_t, 'cm');

%%
load Results_sources_2nd_Day\4-8Hz\sources_Pre_RH_Im1_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta_5new.mat;
load Results_sources_2nd_Day\4-8Hz\sources_Pre_RH_Im2_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta_5new.mat;
load Results_sources_2nd_Day\4-8Hz\sources_Bgr_before_RH_Im1_foi=6Hz_df=2Hz_eLoreta_5new.mat;

subj_set=[1:5];
for subj=subj_set
    source_Norm{subj}=sourcepre(subj);
    source_Norm{subj}.avg.pow=(source_Norm{subj}.avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
    
    source_Norm2{subj}=sourcepre2(subj);
    source_Norm2{subj}.avg.pow=(source_Norm2{subj}.avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
end;
load Stat_results_sources\stat_eloreta_D2_Im1_Base_0.5-4.5s_vs_D2_Im2_Base_0.5-4.5s_RestNorm_f0=6Hz_df=2Hz_CA=0.01_negclust.mat;

%%
clear B C zones_info clust_tissue k;

stat.posmask         = (stat.negclusterslabelmat<2)&(stat.negclusterslabelmat>0);

k=find(stat.mask==1);
  
cfg = [];
cfg.interpmethod='nearest';
cfg.parameter='all';
[int_atlas] = ft_sourceinterpolate(cfg, atlas_t, stat);

for i=1:length(k)
    clust_tissue(i,1)=int_atlas.tissue(k(i));
    clust_tissue(i,2)=k(i);
end;
B=sortrows(clust_tissue,1);

C = unique(B(:,1));

nz=0;
for z=1:length(C)
    if C(z,1)>0
        nz=nz+1;
        zones_info{nz,1}=C(z,1);
        zones_info{nz,2}=atlas.tissuelabel{1,C(z,1)};    
    end;
end;

%%
clear idx;

idx=68;

clustmask         = (B(:,1)==idx); 
clustmask=sum(clustmask, 2); 
k=find(clustmask==1);

cl_size=size(k);
cl_size=cl_size(1);

for subj=subj_set
    p1=0;
    p2=0;
    for iclust=1:cl_size
        source_num=B(k(iclust),2);
        p1=p1+source_Norm{subj}.avg.pow(source_num);
        p2=p2+source_Norm2{subj}.avg.pow(source_num);
    end;
    pow_avg(subj,1)=p1/cl_size;
    pow_avg(subj,2)=p2/cl_size;
end;

save('ClustPow_zone\pow_avg_D2_base1_base2_f0=6Hz_df=2Hz_eloreta_Norm_to_Im1_Bgr_PrecuneusR_in_clust_base1_vs_base2_CA=0.01_5new.txt', 'pow_avg', '-ascii');

clear clustmask;
clear clust_tissue;
clear k;
clear idx;

%%
%pow_avg=load('ClustPow_full_head\pow_avg_D2_Im1_base_post_f0=6Hz_df=2Hz_eloreta_CA=0.0025_Norm_to_Rest.txt', '-ascii');
for n=1:2
    pow_grAvg(n)=0;
    for subj=subj_set
        pow_grAvg(n)=pow_grAvg(n)+pow_avg(subj,n);
    end;
    pow_grAvg(n)=pow_grAvg(n)/5;
end;

figure; plot(pow_grAvg)
