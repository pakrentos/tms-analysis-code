%%
global ft_default
ft_default.trackusage = 'no';
ft_defaults;

% load Results_sources_2nd_Day\4-8Hz\sources_Pre_RH_Im1_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta_RSNorm.mat;
% load Results_sources_2nd_Day\4-8Hz\sources_Pre_RH_Im2_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta_RSNorm.mat;

load Results_sources_2nd_Day\14-30Hz\sources_Pre_RH_Im1_foi=22Hz_df=8Hz_0.5-4.5s_eLoreta_RSNorm.mat;
load Results_sources_2nd_Day\14-30Hz\sources_Post_RH_Im1_foi=22Hz_df=8Hz_6-8s_eLoreta_RSNorm.mat;
pre1 = sourcepre_Norm;
post1 = sourcepost_Norm;

load Results_sources_sham_Day\14-30Hz\sources_Pre_RH_Im1_foi=22Hz_df=8Hz_0.5-4.5s_eLoreta_RSNorm.mat;
load Results_sources_sham_Day\14-30Hz\sources_Post_RH_Im1_foi=22Hz_df=8Hz_6-8s_eLoreta_RSNorm.mat;
pre1(21:35) = sourcepre_Norm;
post1(21:35) = sourcepost_Norm;
 
sourcepre_Norm = pre1;
sourcepost_Norm = post1;

% load Results_sources_2nd_Day\4-8Hz\sources_Pre_RH_Im1_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta.mat;
% load Results_sources_2nd_Day\4-8Hz\sources_Pre_RH_Im2_foi=6Hz_df=2Hz_0.5-4.5s_eLoreta.mat;
% load Results_sources_2nd_Day\4-8Hz\sources_Bgr_before_RH_Im1_foi=6Hz_df=2Hz_eLoreta.mat;

% load Results_sources_2nd_Day\12-14Hz\sources_Bgr_before_RH_Im1_foi=13Hz_df=1Hz_eLoreta.mat;
% sourcebgr1 = sourcebgr;
% clear sourcebgr
% load Results_sources_2nd_Day\12-14Hz\sources_Bgr_before_RH_Im2_foi=13Hz_df=1Hz_eLoreta.mat;
% load Results_sources_2nd_Day\12-14Hz\sources_Bgr_before_RH_Re_foi=13Hz_df=1Hz_eLoreta.mat;


clear C1 C2;
subj_set=[1:9 11:15];
nsubj=0;
for subj=subj_set
    nsubj=nsubj+1;
     C1{nsubj} = sourcepre_Norm(subj);
     C2{nsubj} = sourcepost_Norm(subj);
%     C1{nsubj}=sourcepre(subj);
%     C1{nsubj}.avg.pow=(C1{nsubj}.avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
% 
%     C2{nsubj}=sourcepre2(subj);
%     C2{nsubj}.avg.pow=(C2{nsubj}.avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
% 
%     C1{nsubj}=sourcebgr1(subj);
%     C1{nsubj}.avg.pow=(C1{nsubj}.avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
% 
%     C2{nsubj}=sourcebgr2(subj);
%     C2{nsubj}.avg.pow=(C2{nsubj}.avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
     
    C1{nsubj}.avg = rmfield(C1{nsubj}.avg,'ori');
    C1{nsubj}.avg = rmfield(C1{nsubj}.avg,'mom');
    C2{nsubj}.avg = rmfield(C2{nsubj}.avg,'ori');
    C2{nsubj}.avg = rmfield(C2{nsubj}.avg,'mom');
    C1{nsubj} = rmfield(C1{nsubj},'time');
    C2{nsubj} = rmfield(C2{nsubj},'time');
end;

% for subj=1:15
% %     Pre{subj} = rmfield(Pre{subj},'trialinfo');
% %     Post{subj} = rmfield(Post{subj},'trialinfo');
%     Pre{subj} = rmfield(Pre{subj},'freq');
%     Post{subj} = rmfield(Post{subj},'freq');
% end;

%%

% run statistics over subjects %
cfg=[];
cfg.dim         = C1{1}.dim;
cfg.method      = 'montecarlo';
%cfg.statistic   = 'ft_statfun_actvsblT';
cfg.statistic   = 'ft_statfun_depsamplesT';
%cfg.statistic   = 'ft_statfun_indepsamplesT';
cfg.parameter   = 'pow';
cfg.correctm    = 'cluster';
%cfg.correctm    = 'fdr';
%cfg.numrandomization = 'all';
cfg.numrandomization = 16000;
cfg.alpha       = 0.025; 
cfg.clusteralpha       = 0.02;
cfg.tail        = 0;
cfg.clustertail        = 0;

design  = zeros(2,2*nsubj);
design(1,1:nsubj) = 1;
design(1,nsubj+1:2*nsubj) = 2;
design(2,1:nsubj) = [1:nsubj];
design(2,nsubj+1:2*nsubj) = [1:nsubj];

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;

stat = ft_sourcestatistics(cfg, C1{:}, C2{:});

save('Stat_results_sources\Beta\stat_eloreta_D2_Im1_Base_0.5-4.5s_vs_D2_Im1_Post_6-8s_Im1_RSNorm_f0=22Hz_df=8Hz_CA=0.02_posclust_34subj.mat', 'stat');

%%
mri = ft_read_mri('Subject01.mri');

mri_resliced = ft_volumereslice([], mri);
mri_resliced = ft_convert_units(mri_resliced,'cm');

% % read the atlas
atlas = ft_read_atlas('ROI_MNI_V4.nii');
atlas_t= ft_convert_coordsys(atlas, 'ctf');
atlas_t = ft_convert_units(atlas_t, 'cm');

%%
statplotAll=statplot1;
statplotAll.posmask=statplot1.posmask+statplot2.posmask;
statplotAll.stat=statplot1.stat.*statplot1.posmask+statplot2.stat.*statplot2.posmask;

%%
stat.posmask         = (stat.posclusterslabelmat<2)&(stat.posclusterslabelmat>0);
%stat.posmask         = (stat.negclusterslabelmat==1);
%k=find(stat.posmask==1);

%stat.posmask         = (stat.posclusterslabelmat<4)&(stat.posclusterslabelmat>0);
% stat.posmask         = (stat.posclusterslabelmat==1);
% k=find(statplot.posmask==1);

% cfg = [];
% mri_norm = ft_volumenormalise(cfg, mri);
% 
% mri_rn = ft_volumereslice([], mri_norm);
% mri_rn = ft_convert_units(mri_rn,'cm');

stat.mask2         = (stat.stat>0.75*max(stat.stat));
%stat.mask2         = (stat.stat<0.7*min(stat.stat));
stat.mask = stat.mask2.*stat.posmask;

cfg = [];
cfg.parameter = 'all';
%cfg.downsample = 2;
statplot = ft_sourceinterpolate(cfg, stat, mri_resliced);

cfg = [];
cfg.atlas = atlas_t;
%cfg.roi=atlas_t.tissuelabel(idx);
cfg.method        = 'slice';
cfg.slicerange = [120 200]; % Low Alpha
%cfg.slicerange = [160 200]; % Theta
cfg.nslices = 12;
%cfg.colorbar='no';
%cfg.axis  ='off';
%cfg.funparameter  = 'posclusterslabelmat';
cfg.funparameter  = 'stat';
%cfg.funcolorlim = [-5 0];
cfg.funcolorlim ='zeromax';
cfg.maskparameter = 'posmask';
ft_sourceplot(cfg, statplot);
cfg.location = 'max';
%cfg.locationcoordinates = 'voxel';
%cfg.location = [163 169 174];  %%% Zone of TMS
%cfg.location = [178 174 139];  %%% Zone of interaction cluster*zone
%  cfg.locationcoordinates = 'head';  
%  cfg.location = [-3.1 -0.6 10.2];  %%% Precuneus R point
 % cfg.location = [7.9 4.6 5.1];  %%% TMS
cfg.method        = 'ortho';
ft_sourceplot(cfg, statplot);
