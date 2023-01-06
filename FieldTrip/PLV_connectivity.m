global ft_default
ft_default.trackusage = 'no';
ft_defaults;

f1N=[4 8 12 10 14];
f2N=[8 12 14 14 30];

indx1 = 7312;
indx2 = 22163;

load vol_openmeeg_new;

load elec_cor_new_TMS;

load leadfield_new_TMS;

load Data_In\DAY2_TMS;

subj_set=[1:9 11:15];

for Nfr=1:5
    
    f1=f1N(Nfr);
    f2=f2N(Nfr);
    
    f0=round((f1+f2)/2);
    df=round((f2-f1)/2);
    
    for subj=subj_set
        
        clear data_in data;
        
        data_in(1) = res{1, subj}.right_im1;
        data_in(2) = res{1, subj}.right_im2;
        data_in(3) = res_bgr{1, subj}.right_real;
        data_in(4) = res_bgr{1, subj}.right_im1;
        data_in(5) = res_bgr{1, subj}.right_im2;
        
        for i=1:5
            data_in(i).hdr.elec=elec_cor;
            data_in(i).hdr.nTrials=size(data_in(i).trial,2);
            data_in(i).label=elec_cor.label;
            data_in(i).hdr.label=elec_cor.label;
            data_in(i).hdr.nChans=32;
            data_in(i).hdr.chanunit=elec_cor.chanunit;
            data_in(i).hdr.chantype=elec_cor.chantype;
            data_in(i).elec=elec_cor;
            
            cfg = [];
            cfg.reref               = 'yes';
            cfg.refchannel          = 'all';
            cfg.refmethod = 'avg';
            cfg.demean = 'yes';
            cfg.bpfilter = 'yes';
            cfg.bpfreq        = [f0-df f0+df];
            data_in(i) = ft_preprocessing(cfg, data_in(i));
            
            cfg = [];
            if (i==1)
                cfg.toilim = [0.5 4.5];
                data(i) = ft_redefinetrial(cfg, data_in(i));
                
                cfg.toilim = [5 5.5];
                data(i+5) = ft_redefinetrial(cfg, data_in(i));
                
                cfg.toilim = [6 8];
                data(i+6) = ft_redefinetrial(cfg, data_in(i));
            elseif (i==2)
                cfg.toilim = [0.5 4.5];
                data(i) = ft_redefinetrial(cfg, data_in(i));
                
                cfg.toilim = [5 5.5];
                data(i+6) = ft_redefinetrial(cfg, data_in(i));
                
                cfg.toilim = [6 8];
                data(i+7) = ft_redefinetrial(cfg, data_in(i));
            else
                data(i) = data_in(i);
            end;
            
        end;
        
        clear tlock sources;
        
        cfg = [];
        cfg.covariance='yes';
        for i=1:5
            tlock(i) = ft_timelockanalysis(cfg, data(i));
        end;
        
        cfg=[];
        cfg.method      = 'eloreta';
        cfg.sourcemodel = leadfield;
        cfg.headmodel   = vol;
        cfg.eloreta.fixedori      = 'yes';
        cfg.eloreta.keepfilter  = 'yes';
        cfg.channel = 'all';
        cfg.elec = elec_cor;
        for i=1:5
            sources(i) = ft_sourceanalysis(cfg, tlock(i));
            filt1 = sources(i).avg.filter{1,indx1};
            filt2 = sources(i).avg.filter{1,indx2};
            sourcedata1 = [];
            sourcedata1.label = {'x', 'y', 'z'};
            sourcedata1.time = data(i).time;
            sourcedata2 = [];
            sourcedata2.label = {'x', 'y', 'z'};
            sourcedata2.time = data(i).time;
            for tr=1:length(data(i).trial)
                sourcedata1.trial{tr} = filt1 * data(i).trial{tr}(:,:);
                sourcedata2.trial{tr} = filt2 * data(i).trial{tr}(:,:);
            end;
            
            clear timeseries1 timeseries2;
            
            timeseries1 = cat(2, sourcedata1.trial{:});
            timeseries2 = cat(2, sourcedata2.trial{:});
            
            [u1, s1, v1] = svd(timeseries1, 'econ');
            [u2, s2, v2] = svd(timeseries2, 'econ');
            
            vchannel = [];
            vchannel.label = {'DLPFC_L', 'Precuneus_R'};
            vchannel.time = data(i).time;
            for tr=1:length(data(i).trial)
                vchannel.trial{tr}(1,:) = u1(:,1)' * filt1 * data(i).trial{tr}(:,:);
                vchannel.trial{tr}(2,:) = u2(:,1)' * filt2 * data(i).trial{tr}(:,:);
            end;
            vchan_all{subj,Nfr}(i) = vchannel;
        end;
    end;
end;

for Nfr=1:5
    
    f1=f1N(Nfr);
    f2=f2N(Nfr);
    
    f0=round((f1+f2)/2);
    df=round((f2-f1)/2);
    
    for subj=subj_set
        for i=1:5
            clear freq;
            cfg            = [];
            cfg.output     = 'fourier';
            cfg.method     = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.foi     = f0;
            cfg.tapsmofrq  = df;
            cfg.keeptrials = 'yes';
            cfg.channel    = {'DLPFC_L', 'Precuneus_R'};
            freq = ft_freqanalysis(cfg, vchan_all{subj,Nfr}(i));
            
            cfg = [];
            cfg.method = 'plv';
            conn{subj,Nfr}(i) = ft_connectivityanalysis(cfg, freq);
        end;
    end;
end;

fout=strcat('PLV_EEG_eLoreta_2_virtual_channels.mat');
save(fout, 'conn');