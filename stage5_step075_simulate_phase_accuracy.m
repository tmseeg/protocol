%% Phase accuracy stimulation of non-stimulated data
% Author: Christoph Zrenner

% set paths to relevant toolboxes
addpath(fullfile('toolboxes', 'phastimate'))
addpath(fullfile('toolboxes', 'fieldtrip-20250123'))
ft_defaults

% load data
load(fullfile('data_epoched', 'OPENLOOP_ftdata.mat'));

% remove bad trials and channels manually
data_clean = ft_rejectvisual([], ftdata_eeg_prestim);

cfg = [];
cfg.channel = 'all';
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';
data_clean = ft_preprocessing(cfg, data_clean);

% define and apply spatial filter to be used

montage.labelold = {'C3', 'FC1', 'FC5', 'CP1', 'CP5'};
montage.tra      = [1.00, -0.25, -0.25, -0.25, -0.25];
montage.labelnew = {'C3 SL'};

data_spf = ft_apply_montage(data_clean, montage);

% Causal and acausal estimates of phase at center of epoch

assert(size(data_spf.trial{1}, 1)==1, 'single channel only')

data = cat(1, data_spf.trial{:})'; % [sample x trial]
assert(data_spf.fsample==1000, 'sample rate of 1kHz expected')
D_phastimate = designfilt('bandpassfir', 'FilterOrder', 130, 'CutoffFrequency1', 8, 'CutoffFrequency2', 13, 'SampleRate', data_spf.fsample);

TRANSPORT_DELAY_SAMPLES = 3;

estphase_ix = round(size(data, 1)/2);

causal_data_indices = (-500:-1) + estphase_ix - TRANSPORT_DELAY_SAMPLES;

estphase_causal = phastimate(data(causal_data_indices,:), D_phastimate, 65, 30, 128, TRANSPORT_DELAY_SAMPLES);

D_acausal = designfilt('bandpassiir', 'FilterOrder', 6, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 13, 'SampleRate', data_spf.fsample);

estphase_acausal = hilbert(filtfilt(D_acausal, data));
estphase_acausal = angle(estphase_acausal(estphase_ix,:));

% Visualize distribution of phase using both estimates and phase error

fig = figure(Visible='on');

ax = subplot(2,2,1,polaraxes);
polarhistogram(ax, estphase_causal, BinWidth=pi/12)
ax.ThetaZeroLocation = "top";
title(ax, 'Causal estimates')

ax = subplot(2,2,2,polaraxes);
polarhistogram(ax, estphase_acausal, BinWidth=pi/12)
ax.ThetaZeroLocation = "top";
title(ax, 'Acausal estimates')

ax = subplot(2,2,[3 4],polaraxes);
polarhistogram(ax, estphase_causal-estphase_acausal, BinWidth=pi/24)
ax.ThetaZeroLocation = "top";
title(ax, 'Difference')
