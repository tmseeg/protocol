%Author: Paolo Belardinelli

%% Toolboxes

TOOLBOXPATH = 'toolboxes';
addpath(genpath(fullfile(TOOLBOXPATH, 'hbf_lc_p_v20201022')))
addpath(fullfile(TOOLBOXPATH, '_SOURCE/metafunctions'))
addpath(fullfile(TOOLBOXPATH, 'phastimate'))

%%
to_eliminate={'NAS_pre', 'LPA_pre', 'RPA_pre', 'REF', 'GND','NAS_post', 'LPA_post', 'RPA_post'}
[A ind]= ismember(to_eliminate,elec_realigned.label)

elec_realigned.label(ind) = []
elec_realigned.chanpos(ind,:) = [];
elec_realigned.elecpos(ind,:) = [];
elec_realigned.chantype(ind) = [];
elec_realigned.chanunit(ind) = [];

%%

bmeshesN =struct2cell(headmodel.bmeshes)
bmeshesN=squeeze(bmeshesN);

BmeshesNN{1}.p=bmeshesN(1,1);
BmeshesNN{2}.p=bmeshesN(1,2);
BmeshesNN{3}.p=bmeshesN(1,3);


BmeshesNN{1}.e=bmeshesN(2,1);
BmeshesNN{2}.e=bmeshesN(2,2);
BmeshesNN{3}.e=bmeshesN(2,3);


BmeshesNN{1}.p= cell2mat(BmeshesNN{1}.p)
BmeshesNN{2}.p= cell2mat(BmeshesNN{2}.p)
BmeshesNN{3}.p= cell2mat(BmeshesNN{3}.p)

BmeshesNN{1}.e= cell2mat(BmeshesNN{1}.e)
BmeshesNN{2}.e= cell2mat(BmeshesNN{2}.e)
BmeshesNN{3}.e= cell2mat(BmeshesNN{3}.e)

%% Calculate leadfield
L=mf_calculate_leadfield(BmeshesNN, headmodel.smesh , elec_realigned.label, elec_realigned);

%%getting the LAMBDA and the COV matrix
lambda_factor=10;
cfg = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = [-1 -0.05];
cfg.keeptrials = 'no';
TEP= ft_timelockanalysis(cfg,ftdata_eeg_prestim);
Cov = TEP.cov;
lambda=lambda_factor*max(eig(Cov));
invCy = pinv(Cov + lambda * eye(size(Cov)));

% Interactively Select M1 Handknob

[M1_center, M1_dipoles, M1_indexes]=mf_select_dipole(headmodel.smesh, BmeshesNN, 2)
lf1 = L(:,find(M1_indexes));
filt = pinv(lf1' * invCy * lf1) * lf1' * invCy;  
source_M1=cellfun(@(x) filt * x, ftdata_eeg_prestim.trial, 'un', 0)';

%Interactively Select S1 Handknob

[S1_center, S1_dipoles, S1_indexes]=mf_select_dipole(headmodel.smesh, BmeshesNN, 2)
lf1 = L(:,find(S1_indexes));
filt = pinv(lf1' * invCy * lf1) * lf1' * invCy;  
source_S1=cellfun(@(x) filt * x, ftdata_eeg_prestim.trial, 'un', 0)';
 
%% Calculate Hjorth sensor Signal around C3

hjorth_spot=ftdata_eeg_prestim;
hjorth_spot.trial={};
hjorth_spot.label={'Hjorth'};
hjort_vector=zeros(length(ftdata_eeg_prestim.label),1);
hjort_vector(find(strcmp(ftdata_eeg_prestim.label, 'C3')))=1;
hjort_vector(find(strcmp(ftdata_eeg_prestim.label, 'FC1')))=-1/4;
hjort_vector(find(strcmp(ftdata_eeg_prestim.label, 'FC5')))=-1/4;
hjort_vector(find(strcmp(ftdata_eeg_prestim.label, 'CP1')))=-1/4;
hjort_vector(find(strcmp(ftdata_eeg_prestim.label, 'CP5')))=-1/4;

for i=1:length(ftdata_eeg_prestim.trial)
    hjorth_spot.trial{i}=hjort_vector'*ftdata_eeg_prestim.trial{i};
end

M1_spot=ftdata_eeg_prestim;
M1_spot.trial={};
M1_spot.label={'Source'};
M1_spot.trial=source_M1';

S1_spot=ftdata_eeg_prestim;
S1_spot.trial={};
S1_spot.label={'Source'};
S1_spot.trial=source_S1';

%%
for i=1:size(S1_spot.trial, 2)
S1_spot.trial{i}=mean(S1_spot.trial{i}, 1)
end

for i=1:size(M1_spot.trial, 2)
M1_spot.trial{i}=mean(M1_spot.trial{i}, 1)
end

%% Calculate invidual spectra for Hjorth, M1, S1. Plot them

cfg=[];
cfg.latency=[-1 -0.005];
hjorth_prestim=ft_selectdata(cfg, hjorth_spot);
M1_prestim=ft_selectdata(cfg, M1_spot);
S1_prestim=ft_selectdata(cfg, S1_spot);

cfg = [];
cfg.output = 'pow';
cfg.channel = 'all';
cfg.method = 'mtmfft';
cfg.foilim  = [2 45];
cfg.tapsmofrq = 2;
cfg.taper = 'hanning';
cfg.keeptrials = 'no';
hjorth_spect = ft_freqanalysis(cfg, hjorth_prestim);
M1_spect = ft_freqanalysis(cfg, M1_prestim);
S1_spect = ft_freqanalysis(cfg, S1_prestim);



lw=4;
plot(hjorth_spect.powspctrm, 'LineWidth', lw)
xlim([2 45])
set(gca,'FontSize',18,'FontWeight','bold','TickLength',[0 0])
ylabel('Absolute power (uV^2) ','FontSize',18,'FontWeight','bold')
xlabel('Frequency (Hz) ','FontSize',18,'FontWeight','bold')
set(gca,'linewidth',2)

figure(2)
plot(M1_spect.powspctrm, 'LineWidth', lw)
xlim([2 45])
set(gca,'FontSize',18,'FontWeight','bold','TickLength',[0 0])
ylabel('a.u.','FontSize',18,'FontWeight','bold')
xlabel('Frequency (Hz) ','FontSize',18,'FontWeight','bold')
set(gca,'linewidth',2)
figure(3)
plot(S1_spect.powspctrm, 'LineWidth', lw)
xlim([2 45])
set(gca,'FontSize',18,'FontWeight','bold','TickLength',[0 0])
ylabel('a.u.','FontSize',18,'FontWeight','bold')
xlabel('Frequency (Hz) ','FontSize',18,'FontWeight','bold')
set(gca,'linewidth',2)


