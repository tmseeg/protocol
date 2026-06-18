%% Post-hoc validation of phase accuracy
% Author: Christoph Zrenner

% set paths to relevant toolboxes
addpath(fullfile('toolboxes', 'fieldtrip-20250123'))
ft_defaults

% load data
...


% create epochs
...


% downsample to 1kHz
...


% demean epoch
...

% remove bad trials and channels manually
data_clean = ft_rejectvisual([], ftdata_eeg_prestim);

% define and apply spatial filter to be used
montage.labelold = {'C3', 'FC1', 'FC5', 'CP1', 'CP5'};
montage.tra      = [1.00, -0.25, -0.25, -0.25, -0.25];
montage.labelnew = {'C3 SL'};
data_spf = ft_apply_montage(data_clean, montage);

%% Step 1: Spectrum

cfg = [];
cfg.foilim        = [3 60];
cfg.pad     = 'nextpow2';
cfg.output  = 'pow';
cfg.channel = 'all';
cfg.method  = 'mtmfft';
cfg.taper   = 'hanning';
spect       = ft_freqanalysis(cfg, data_spf);

figure; ax = axes;
plot(ax, spect.freq, pow2db(spect.powspctrm));
ax.XScale = "log";

cfg               = [];
cfg.channel       = {'C3'};
cfg.foilim        = [3 60];
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
cfg.output        = 'fractal';
fractal = ft_freqanalysis(cfg, data_spf);
cfg.output        = 'original';
original = ft_freqanalysis(cfg, data_spf);

% remove the aperiodic part from the full spectrum
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x1/x2';
oscillatory = ft_math(cfg, original, fractal);

figure;

ax1 = subplot(1,2,1); hold(ax1, "on");
plot(original.freq, pow2db(original.powspctrm),'k');
plot(fractal.freq, pow2db(fractal.powspctrm));
legend({'original','aperiodic'},'location','southwest');
ax1.XScale = "log";
xlabel('Freq (Hz)'); ylabel('PSD (dB)');

ax2 = subplot(1,2,2); hold(ax2, "on");
plot(oscillatory.freq, pow2db(oscillatory.powspctrm));
ax2.XScale = "log";
xlabel('Freq (Hz)'); ylabel('SNR (dB)');


% Post-hoc estimate of phase at center of epoch (non-stimulated data, no TMS artifact)

assert(size(data_spf.trial{1}, 1)==1, 'single channel only')

data = cat(1, data_spf.trial{:})'; % [sample x trial]
estphase_ix = round(size(data, 1)/2);

assert(data_spf.fsample==1000, 'sample rate of 1kHz expected')

% NOTE: to avoid reproducing filter-related biases, it's recommended to use a different filter design than what is used in the causal real-time estimate, where an FIR filter is used
D = designfilt('bandpassiir', 'FilterOrder', 6, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 13, 'SampleRate', data_spf.fsample);

estphase_acausal = hilbert(filtfilt(D, data_spf));
estphase_acausal = angle(estphase_acausal(estphase_ix,:));

% Visualize distribution of phase using both estimates and phase error

fig = figure(Visible='on');

ax = subplot(2,2,2,polaraxes);
polarhistogram(ax, estphase_acausal, BinWidth=pi/12)
ax.ThetaZeroLocation = "top";
title(ax, 'Actual phase at marker')

