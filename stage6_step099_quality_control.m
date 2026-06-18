%% Quality Control
% Author: Christoph Zrenner

% set paths to relevant toolboxes
addpath(fullfile('toolboxes', 'eeglab-2024.2.1')) % with the neurone plug-in installed, std_maketrialinfo.m is patched to work with eeglab2fieldtrip, line 164 replaced with if ~isempty(STUDY) & isfield(STUDY.datasetinfo, 'trialinfo')
eeglab; % launch eeglab

% load epoched data around TMS pulse
load(fullfile('data_epoched', 'REFTEP_001_eeglab.mat'));

% Note that REFTEP_001 data is TMS data collected in an open-loop fashion

EVENT_TYPES = {'A - Out'};

assert(size(EEG.data, 1) == 5, 'this script expects 5-channel EEG-data with the first channel in the centre')

%% Apply Hjorth-style Surfac Laplacian Spatial Filter

EEG_spf = EEG;
EEG_spf.data = pagemtimes([1 -0.25 -0.25 -0.25 -0.25], EEG.data);
eeg_checkset(EEG_spf);
EEG_spf.chanlocs(1).labels = {'Hjorth SL'};

% Better code:
% [chan_ismember, chan_ix] = ismember({'C3', 'FC1', 'FC5', 'CP1', 'CP5'}, {EEG.chanlocs.labels})
% assert(all(chan_ismember), 'EEG data must contain all channels used in the spatial filter')
% EEG.data = cat(1, EEG.data, pagemtimes([1 -0.25 -0.25 -0.25 -0.25], EEG.data(chan_ix, :, :)))
% EEG.chanlocs(end+1).labels = 'Hjorth SL'; EEG.nbchan = EEG.nbchan+1;
% eeg_checkset(EEG);

%% Extract and downsample pre-stimulus period

EEG_prestimulus = EEG;

EEG_prestimulus = pop_select(EEG_prestimulus, 'time', [-1.410 -0.010]);
EEG_prestimulus = pop_resample(EEG_prestimulus, 1000); 
EEG_prestimulus = pop_rmbase(EEG_prestimulus, [], []);


%% Spectrum from epochs

SPECT_EPOCHLENGTH = 1.4; % seconds
SPECT_EPOCHOVERLAP = 0.5 * SPECT_EPOCHLENGTH; % seconds
SPECT_DPSS_NW = 4/2;

FREQBIN_MIN = 2; % Hz
FREQBIN_MAX = 70; % Hz
FREQBIN_LOGSPACENUM = 100; % duplicate bins may be removed

epochs = [];
epochs.trial = double(permute(EEG_prestimulus.data, [2 1 3]));
epochs.dimord = 'time_chan_rpt'; % fieldtrip uses 'rpt_chan_time', eeglab uses 'chan_time_rpt', our analysis code expects 'time_chan_rpt'
epochs.time = EEG_prestimulus.times/1000; % convert from ms to s
epochs.label = {EEG_prestimulus.chanlocs.labels}';


%% Perform Spectral Analysis

Fs = round(1/mean(diff(epochs.time)));
assert(Fs == 1000, 'expecting sample rate of 1 kHz');

irasa_factors = 1.1:0.1:2.9; % see also https://doi.org/10.1101/203786
irasa_factors(irasa_factors == 2) = []; %exclude factor 2.0

nfft = max(irasa_factors) * size(epochs.trial, 1); % minimum padding for IRASA resampling
nfft = 2^nextpow2(nfft); % radix-2 for more efficient fft

nw = SPECT_DPSS_NW; % 4/2, 5/2, 6/2, 7/2; % time-halfbandwidth product
ntp = 2*(nw)-1; % number of tapers

freqboi_cf = Fs*(0:(nfft/2))/nfft; % center frequency for each bin, including DC bin  
%freqboi_ix = find((freqboi_cf >= FREQBIN_MIN) & (freqboi_cf <= FREQBIN_MAX)); % linear
[~, freqboi_ix] = min(abs(logspace(log10(FREQBIN_MIN), log10(FREQBIN_MAX), FREQBIN_LOGSPACENUM) - freqboi_cf')); % log
freqboi_ix = unique(freqboi_ix);

freqboi_cf = freqboi_cf(freqboi_ix);

tprs = dpss(size(epochs.trial, 1), nw, ntp); % multi-taper
%tprs(:, end) = []; % drop last taper

%tpr = hanning(size(epochs.trial, 1)); % window function

% tpr = tpr ./ norm(tpr, 'fro') % does this mess up the density correction? % Note: apply same procedure in IRASA

fprintf('\nSpectral estimation parameters: nfft=%i, ntp=%i ...', nfft, ntp)

spect = [];
spect.fourierspctrm = fft(epochs.trial .* permute(tprs, [1 3 4 2]), nfft);
spect.fourierspctrm = spect.fourierspctrm(freqboi_ix,:,:,:) .* sqrt(2 ./ nfft); % density and 1-sided spectrum correction
spect.fourierspctrm = permute(spect.fourierspctrm, [2 1 3 4]);
spect.freq = freqboi_cf;
spect.label = epochs.label;
spect.dimord = 'chan_freq_rpt_tpr'; % [channel, bin, trial, taper]

fprintf(' done.\n')


%% Perform by-trial IRASA

assert(exist('epochs') == 1); % epochs.trial should be [time, channel, trial]
assert(exist('irasa_factors') == 1); % irasa resampling factors

assert((max(irasa_factors) * size(epochs.trial, 1)) < nfft, 'insufficient zero padding for all resampling factors')

fprintf('\nPerforming sensor level IRASA with %i factors ..', numel(irasa_factors)), tic;

% 5-d matrix with [channel, bin, trial, taper, resampling factor]
ucom_ih = complex(zeros(size(epochs.trial, 2), length(freqboi_ix), size(epochs.trial, 3), ntp, length(irasa_factors), 'single'));
dcom_ih = complex(zeros(size(ucom_ih), 'single'));

for ih = 1:length(irasa_factors) % loop across resampling factors
    fprintf('.');

    % resample sensor timecourses and apply tapers
    [n, d] = rat(irasa_factors(ih)); % n > d

    udat = resample(reshape(epochs.trial, size(epochs.trial, 1), []), n, d); % upsample
    udat = reshape(udat, [], size(epochs.trial, 2), size(epochs.trial, 3));

    %utpr = hanning(size(udat, 1)); %utpr = utpr ./ norm(utpr, 'fro');
    %udat = udat .* utpr; % apply taper

    utpr = dpss(size(udat, 1), nw, ntp); % compute tapers
    udat = udat .* permute(utpr, [1 3 4 2]); % apply taper

    ddat = resample(reshape(epochs.trial, size(epochs.trial, 1), []), d, n); % downsample
    ddat = reshape(ddat, [], size(epochs.trial, 2), size(epochs.trial, 3));

    dtpr = dpss(size(ddat, 1), nw, ntp); % compute tapers
    ddat = ddat .* permute(dtpr, [1 3 4 2]); % apply taper

    % fft of upsampled data [bin x trial x channel x taper]
    ucom = fft(udat, nfft);
    ucom = ucom(freqboi_ix,:,:,:);
    ucom = ucom .* sqrt(2 ./ nfft);

    % fft of downsampled data
    dcom = fft(ddat, nfft);
    dcom = dcom(freqboi_ix,:,:,:);
    dcom = dcom .* sqrt(2 ./ nfft);

    ucom_ih(:,:,:,:,ih) = permute(single(ucom), [2 1 3 4]);
    dcom_ih(:,:,:,:,ih) = permute(single(dcom), [2 1 3 4]);
end
fprintf(' done (took %1.1f s, yielded two %1.1f MB matrices).\n', toc, getfield(whos('ucom_ih'), 'bytes') / 1e6);


%% Apply spatial filter transforms to compute spectra

fprintf('Computing IRASA transforms in frequency domain:\n')
w = [[1 -0.25 -0.25 -0.25 -0.25]; eye(5)];
    
% 4-d matrix with [channel, bin, trial, taper]
fourierspctrm_transformed = pagemtimes(w, spect.fourierspctrm);
original_powspctrm = mean(abs(fourierspctrm_transformed).^2, 4); % mean over tapers

% 5-d matrix with [channel, bin, trial, resampling factor, taper]
ucom_ih_transformed = pagemtimes(w, ucom_ih);
dcom_ih_transformed = pagemtimes(w, dcom_ih);
               
% now take the geometric mean of the resampling factors
% note that after IRASA we have a powspctrm and no longer a complex-valued fourierspctrm
% [bin x channel x trial]
irasa_powspctrm = mean(median(sqrt((abs(ucom_ih_transformed).^2) .* (abs(dcom_ih_transformed).^2)), 4), 5); % (1) power abs()^2,  (2) geometric mean of u and d sqrt(a*b), (3) median over resamlpling factors, (5) mean over tapers

assert(~any(irasa_powspctrm(:) == 0), 'no entries in the 1/f powerspectrum can be zero')

% Compute spectral estimates
psd_original_mean      =    mean(original_powspctrm, 3);
psd_original_prctile25 = prctile(original_powspctrm, 25, 3);        
psd_original_prctile50 = prctile(original_powspctrm, 50, 3);
psd_original_prctile75 = prctile(original_powspctrm, 75, 3);

psd_fractal_mean      =    mean(irasa_powspctrm, 3);
psd_fractal_prctile25 = prctile(irasa_powspctrm, 25, 3);
psd_fractal_prctile50 = prctile(irasa_powspctrm, 50, 3);
psd_fractal_prctile75 = prctile(irasa_powspctrm, 75, 3);        

psd_snr_mean      =    mean(original_powspctrm ./ irasa_powspctrm, 3);
psd_snr_prctile25 = prctile(original_powspctrm ./ irasa_powspctrm, 25, 3);
psd_snr_prctile50 = prctile(original_powspctrm ./ irasa_powspctrm, 50, 3);
psd_snr_prctile75 = prctile(original_powspctrm ./ irasa_powspctrm, 75, 3);


%% Visualize

if ischar(EEG.urevent(1).type)
    ev_times = [EEG.urevent(ismember({EEG.urevent.type}, EVENT_TYPES)).latency]./EEG.srate;
else
    ev_times = [EEG.urevent(ismember([EEG.urevent.type], [EVENT_TYPES{:}])).latency]./EEG.srate;
end
%%
fig = figure('Color', 'white', Name="Quality Control", Units='in', PaperOrientation='portrait', Position=2.*[0 0 8.5 11]); % letter
t = tiledlayout(fig, 'flow');

ax = nexttile(t);
plot(ax, ev_times, '.')
xlabel('Stimulus #')
ylabel('Time (s)')
text(ax, 0, ev_times(end), sprintf(' %d trials, %d seconds', numel(ev_times), round(ev_times(end))), HorizontalAlignment="left", VerticalAlignment="top")

ax = nexttile(t);

histogram(ax, log10(diff(ev_times)), BinWidth=0.1)
ax.XLim(1) = 0;
ax.XTick = log10([1 2 5 10 20])
ax.XTickLabel = {'1', '2', '5', '10', '20'}
text(ax, 0, ax.YLim(2), sprintf(' median=%.1f s\n max=%.1f s', median(diff(ev_times)), max(diff(ev_times))), HorizontalAlignment="left", VerticalAlignment="top")
title('Stimulus Interval Distribution')

ax = nexttile(t);
plot(ax, EEG.times, mean(EEG.data, 3)')
ylim(ax, [-20 20]), ylabel(ax, "µV")
xlim(ax, [-300 300]), xlabel(ax, "ms")
legend(ax, {EEG.chanlocs.labels}, Location="northwest")
title(sprintf('average TEP, original reference (n=%i)', size(EEG.data, 3)))

for trial_ix = floor(linspace(1,size(EEG.data,3),11))
ax = nexttile(t); hold(ax, "on")
plot(ax, EEG.times, EEG.data(:,:,trial_ix))
plot(ax, EEG_spf.times, EEG_spf.data(:,:,trial_ix), 'b', LineWidth=1)
ylim(ax, [-50 50]), ylabel(ax, "µV")
xlim(ax, [-250 250]), xlabel(ax, "ms")
title(sprintf('single trial #%i, original and SL', trial_ix))
end

ax = nexttile(t); hold(ax, "on")
%plot(ax, EEG_spf.times, squeeze(EEG_spf.data(1,:,:)) + ((1:size(EEG_spf.data, 3))-size(EEG_spf.data, 3)/2)*10)
%plot(ax, EEG_spf.times, EEG_spf.data(:,:,trial_ix))
plot(ax, EEG_spf.times, mean(EEG_spf.data, 3)', LineWidth=1.5)
ylim(ax, [-10 10]), ylabel(ax, "µV")
xlim(ax, [-250 50]), xlabel(ax, "ms")
title('Surface Laplacian, average pre-stimulus')

ax = nexttile(t); hold(ax, "on")

timewindow_mask = EEG.times > -0.500 & EEG.times < -0.010;
trials_range = mean(range(EEG_spf.data(:,timewindow_mask,:), 2), 1);
trials_range = trials_range(:);

plot(ax, EEG_spf.times, mean(EEG_spf.data(:,:,trials_range < prctile(trials_range, 80)), 3)', LineWidth=1.5)
ylim(ax, [-10 10]), ylabel(ax, "µV")
xlim(ax, [-250 50]), xlabel(ax, "ms")
title('Surface Laplacian, average pre-stimulus (80% of trials, by range)')

ax = nexttile(t); hold(ax, "on")

timewindow_mask = EEG.times > -0.500 & EEG.times < -0.010;

plot(ax, EEG_spf.times, mean(EEG_spf.data(:,:,1:ceil(size(EEG_spf.data, 3)/3)), 3)', LineWidth=1.5)
ylim(ax, [-10 10]), ylabel(ax, "µV")
xlim(ax, [-250 50]), xlabel(ax, "ms")
title('Surface Laplacian, average pre-stimulus (first 33% of trials)')

ax = nexttile(t); hold(ax, 'on')
%plot(ax, spect.freq, pow2db(mean(abs(spect.fourierspctrm).^2, [3 4])))
plot(ax, spect.freq, pow2db(psd_original_mean))
ax.ColorOrderIndex = 1;
plot(ax, spect.freq, pow2db(psd_fractal_mean), LineStyle=":")
ax.Children(6).LineWidth = 2.5;
ax.Children(12).LineWidth = 2.5;
legend(ax, {'Surface Laplacian' spect.label{:}})
set(ax, XScale='log')
xlim(ax, [4 55])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title('PSD')

ax = nexttile(t);
plot(ax, spect.freq, pow2db(psd_snr_mean))
ax.Children(6).LineWidth = 2.5;
set(ax, XScale='log')
xlim(ax, [4 55])
xlabel('Frequency (Hz)')
ylabel('SNR (dB)')
title('1/f corrected spectrum')

ax = nexttile(t);
plot(ax, EMG.times, mean(EMG.data, 3)')
%ylim(ax, [-1e3 1e3])
ylabel(ax, "µV")
xlim(ax, [-25 75]), xlabel(ax, "ms")
title('Average MEP')