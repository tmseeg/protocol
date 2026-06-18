%% Spectral Estimation using Fieldtrip
% Author: Christoph Zrenner

% set paths to relevant toolboxes
addpath(fullfile('toolboxes', 'fieldtrip-20250123'))

ft_defaults

load(fullfile('data_epoched', 'OPENLOOP_ftdata.mat'));

data_clean = ft_rejectvisual([], ftdata_eeg_prestim);

% remove bad trials and channels manually

% Re-reference to CAR (can be skipped, as choice of reference does not affect subsequent transformation)
cfg = [];
cfg.channel = 'all';
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';
data = ft_preprocessing(cfg, data_clean);

% Option 1: Apply SCD-style Surface Laplacian across all senors
%data = ft_scalpcurrentdensity([], data);

% Option 2: Apply Hjorth-style Surface Laplacian to channel C3 only
data = ft_apply_montage(data, struct('labelold', {{'C3', 'FC1', 'FC5', 'CP1', 'CP5'}}, 'tra', [1.00, -0.25, -0.25, -0.25, -0.25], 'labelnew', {{'C3'}}));


cfg = [];
cfg.foilim        = [3 60];
cfg.pad     = 'nextpow2';
cfg.output  = 'pow';
cfg.channel = 'all';
cfg.method  = 'mtmfft';
cfg.taper   = 'hanning';

spect       = ft_freqanalysis(cfg, data);

figure; ax = axes;
plot(ax, spect.freq, pow2db(spect.powspctrm));
ax.XScale = "log";

cfg               = [];
cfg.channel       = {'C3'};
cfg.foilim        = [3 60];
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
cfg.output        = 'fractal';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'original';
original = ft_freqanalysis(cfg, data);

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

