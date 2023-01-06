%%
global ft_default
ft_default.trackusage = 'no';
ft_defaults;

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
  
    data_in_bgr=res_bgr{1, subj}.right_im1;
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
    
    data_in_bgr2=res_bgr{1, subj}.right_im2;
    for i=1:32
        data_in_bgr2.hdr.chanunit{i,1}='V';
        data_in_bgr2.hdr.chantype{i,1}='eeg';
    end;
    data_in_bgr2.hdr.chanunit(33:34,:)=[];
    data_in_bgr2.hdr.chantype(33:34,:)=[];
    data_in_bgr2.label = elec_cor.label;
    data_in_bgr2.hdr.label = data_in_bgr2.label;
    data_in_bgr2.hdr.nChans=32;
    data_in_bgr2.hdr.elec=elec_cor;
    data_in_bgr2.hdr.nTrials=size(data_in_bgr2.trial,2);
    data_in_bgr2.elec=data_in_bgr2.hdr.elec;
    
    cfg = [];
    cfg.reref               = 'yes';  
    cfg.refchannel          = 'all';    
    cfg.demean = 'yes';
    cfg.bpfilter = 'yes';
    cfg.bpfreq        = [f0-df f0+df];
    data_eeg_bgr_reref          = ft_preprocessing(cfg, data_in_bgr);
    data_eeg_bgr2_reref          = ft_preprocessing(cfg, data_in_bgr2);

    cfg = [];
    cfg.covariance='yes';
    tlock_bgr = ft_timelockanalysis(cfg,data_eeg_bgr_reref);
    tlock_bgr.elec = elec_cor;
    tlock_bgr2 = ft_timelockanalysis(cfg,data_eeg_bgr2_reref);
    tlock_bgr2.elec = elec_cor;
        
    cfg=[];
    cfg.method='eloreta';
    cfg.grid=leadfield;
    cfg.elec = elec_cor;
    cfg.headmodel=vol;
    cfg.channel = 'all';
    cfg.senstype = 'EEG';
    sourcebgr(subj)=ft_sourceanalysis(cfg, tlock_bgr);
    sourcebgr2(subj)=ft_sourceanalysis(cfg, tlock_bgr2);
    
    sourcebgr(subj).avg.mom=[];
    sourcebgr(subj).avg.ori=[];
    sourcebgr2(subj).avg.mom=[];
    sourcebgr2(subj).avg.ori=[];

end;

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Bgr_before_RH_Im1_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_eLoreta.mat');
save(fname, 'sourcebgr');

fname=strcat('Results_sources_2nd_Day\',int2str(f1),'-',int2str(f2),'Hz\sources_Bgr_before_RH_Im2_foi=',int2str(f0),'Hz_df=',int2str(df),'Hz_eLoreta.mat');
save(fname, 'sourcebgr2');

clear data_in_bgr data_eeg_bgr_reref data_in_bgr2 data_eeg_bgr2_reref;
clear tlock_bgr sourcebgr tlock_bgr2 sourcebgr2;

end;
