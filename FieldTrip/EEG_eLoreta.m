%%
global ft_default
ft_default.trackusage = 'no';
ft_defaults;

ft_prepare_leadfield

t1=5;
t2=5.5;

t1_base=0.5;
t2_base=4.5;

load vol_openmeeg_new;

load elec_cor_new_TMS;

load leadfield_new_TMS;


load Data_In\DAY2_TMS;
%load Data_In\DAY2_SHAM;

subj_set=[1:9 11:15];

f1N=[4 8 14 12 10];
f2N=[8 12 30 14 14];

for Nfr=5:5

f1=f1N(Nfr);
f2=f2N(Nfr);

f0=round((f1+f2)/2);
df=round((f2-f1)/2);    
    
for subj=subj_set
    data_in=res{1, subj}.right_im1;
    for i=1:32
        data_in.hdr.chanunit{i,1}='V';
        data_in.hdr.chantype{i,1}='eeg';
    end;
    data_in.hdr.chanunit(33:34,:)=[];
    data_in.hdr.chantype(33:34,:)=[];
    data_in.label = elec_cor.label;
    data_in.hdr.label = data_in.label;
    data_in.hdr.nChans=32;
    data_in.hdr.elec=elec_cor;
    data_in.hdr.nTrials=size(data_in.trial,2);
    data_in.elec=data_in.hdr.elec;
    
    data_in2=res{1, subj}.right_im2;
    for i=1:32
        data_in2.hdr.chanunit{i,1}='V';
        data_in2.hdr.chantype{i,1}='eeg';
    end;
    data_in2.hdr.chanunit(33:34,:)=[];
    data_in2.hdr.chantype(33:34,:)=[];
    data_in2.label = elec_cor.label;
    data_in2.hdr.label = data_in2.label;
    data_in2.hdr.nChans=32;
    data_in2.hdr.elec=elec_cor;
    data_in2.hdr.nTrials=size(data_in2.trial,2);
    data_in2.elec=data_in2.hdr.elec;

    data_in_bgr=res_bgr{1, subj}.right_real;
    for i=1:32
        data_in_bgr.hdr.chanunit{i,1}='V';
        data_in_bgr.hdr.chantype{i,1}='eeg';
    end;
    data_in_bgr.hdr.chanunit(33:34,:)=[];
    data_in_bgr.hdr.chantype(33:34,:)=[];
    data_in_bgr.label = elec_cor.label;
    data_in_bgr.hdr.label = data_in_bgr.label;
    data_in_bgr.hdr.nChans=32;
    data_in_bgr.hdr.elec=elec_cor;
    data_in_bgr.hdr.nTrials=size(data_in_bgr.trial,2);
    data_in_bgr.elec=data_in_bgr.hdr.elec;
    
    cfg = [];
    cfg.reref               = 'yes';   
    cfg.refchannel          = 'all';   
    cfg.refmethod     = 'avg';
    cfg.demean = 'yes';
    cfg.bpfilter = 'yes';
    cfg.bpfreq        = [f0-df f0+df];
    data_eeg_reref          = ft_preprocessing(cfg, data_in);
    data_eeg_reref2          = ft_preprocessing(cfg, data_in2);
    data_eeg_bgr_reref          = ft_preprocessing(cfg, data_in_bgr);
   
    cfg = [];
    cfg.toilim = [t1_base t2_base];
    datapre = ft_redefinetrial(cfg, data_eeg_reref);
    datapre2 = ft_redefinetrial(cfg, data_eeg_reref2);
    
    cfg = [];
    cfg.toilim = [t1 t2];
    datapost = ft_redefinetrial(cfg, data_eeg_reref);
    datapost2 = ft_redefinetrial(cfg, data_eeg_reref2);

    cfg = [];
    cfg.covariance='yes';
    tlock_pre = ft_timelockanalysis(cfg,datapre);
    tlock_pre.elec = elec_cor;
    tlock_post = ft_timelockanalysis(cfg,datapost);
    tlock_post.elec = elec_cor;
    tlock_pre2 = ft_timelockanalysis(cfg,datapre2);
    tlock_pre2.elec = elec_cor;
    tlock_post2 = ft_timelockanalysis(cfg,datapost2);
    tlock_post2.elec = elec_cor;
    tlock_bgr = ft_timelockanalysis(cfg,data_eeg_bgr_reref);
    tlock_bgr.elec = elec_cor;

    cfg=[];
    cfg.method='eloreta';
    cfg.grid=leadfield;
    cfg.elec = elec_cor;
    cfg.headmodel=vol;
    cfg.channel = 'all';
    cfg.senstype = 'EEG';
    sourcepre(subj)=ft_sourceanalysis(cfg, tlock_pre);
    sourcepost(subj)=ft_sourceanalysis(cfg, tlock_post);
    sourcepre2(subj)=ft_sourceanalysis(cfg, tlock_pre2);
    sourcepost2(subj)=ft_sourceanalysis(cfg, tlock_post2);
    sourcebgr(subj)=ft_sourceanalysis(cfg, tlock_bgr);
    
    sourcepre(subj).avg.mom=[];
    sourcepre(subj).avg.ori=[];
    sourcepost(subj).avg.mom=[];
    sourcepost(subj).avg.ori=[];
    sourcepre2(subj).avg.mom=[];
    sourcepre2(subj).avg.ori=[];
    sourcepost2(subj).avg.mom=[];
    sourcepost2(subj).avg.ori=[];
    sourcebgr(subj).avg.mom=[];
    sourcebgr(subj).avg.ori=[];
    
    sourcepre_Norm(subj)=sourcepre(subj);
    sourcepre_Norm(subj).avg.pow=(sourcepre_Norm(subj).avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
    
    sourcepost_Norm(subj)=sourcepost(subj);
    sourcepost_Norm(subj).avg.pow=(sourcepost_Norm(subj).avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
    
    sourcediff(subj)=sourcepost(subj);
    sourcediff(subj).avg.pow=(sourcediff(subj).avg.pow-sourcepre(subj).avg.pow)./sourcepre(subj).avg.pow;
    
    sourcepre_Norm2(subj)=sourcepre2(subj);
    sourcepre_Norm2(subj).avg.pow=(sourcepre_Norm2(subj).avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
    
    sourcepost_Norm2(subj)=sourcepost2(subj);
    sourcepost_Norm2(subj).avg.pow=(sourcepost_Norm2(subj).avg.pow-sourcebgr(subj).avg.pow)./sourcebgr(subj).avg.pow;
    
    sourcediff2(subj)=sourcepost2(subj);
    sourcediff2(subj).avg.pow=(sourcediff2(subj).avg.pow-sourcepre2(subj).avg.pow)./sourcepre2(subj).avg.pow;

end;

fname=strcat('Results_sources_sham_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Pre_RH_Im1_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1_base),'-',num2str(t2_base),'s_eLoreta_RSNorm.mat');
save(fname, 'sourcepre_Norm');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Post_RH_Im1_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1),'-',num2str(t2),'s_eLoreta_RSNorm.mat');
save(fname, 'sourcepost_Norm');
 
fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Pre_RH_Im2_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1_base),'-',num2str(t2_base),'s_eLoreta_RSNorm.mat');
save(fname, 'sourcepre_Norm2');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Post_RH_Im2_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1),'-',num2str(t2),'s_eLoreta_RSNorm.mat');
save(fname, 'sourcepost_Norm2');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Pre_RH_Im1_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1_base),'-',num2str(t2_base),'s_eLoreta.mat');
save(fname, 'sourcepre');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Post_RH_Im1_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1),'-',num2str(t2),'s_eLoreta.mat');
save(fname, 'sourcepost');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Pre_RH_Im2_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1_base),'-',num2str(t2_base),'s_eLoreta.mat');
save(fname, 'sourcepre2');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Post_RH_Im2_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1),'-',num2str(t2),'s_eLoreta.mat');
save(fname, 'sourcepost2');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Bgr_before_RH_Re_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_eLoreta.mat');
save(fname, 'sourcebgr');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Diff_Post-Pre_RH_Im1_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1),'-',num2str(t2),'s_eLoreta.mat');
save(fname, 'sourcediff');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Diff_Post-Pre_RH_Im2_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_',num2str(t1),'-',num2str(t2),'s_eLoreta.mat');
save(fname, 'sourcediff2');    

clear data_in data_in2 data_in_bgr data_eeg_reref data_eeg_reref2 data_eeg_bgr_reref datapre datapre2 datapost datapost2;
clear tlock_pre tlock_pre2 tlock_post tlock_post2 tlock_bgr sourcepre sourcepre2 sourcepost sourcepost2 sourcebgr;
clear sourcepre_Norm sourcepre_Norm2 sourcepost_Norm sourcepost_Norm2 sourcediff sourcediff2;

end;

