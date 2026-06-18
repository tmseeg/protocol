%% Create Epochs

% set paths to relevant toolboxes
addpath(fullfile('toolboxes', 'eeglab-2024.2.1')) % with the neurone plug-in installed, std_maketrialinfo.m is patched to work with eeglab2fieldtrip, line 164 replaced with if ~isempty(STUDY) & isfield(STUDY.datasetinfo, 'trialinfo')
eeglab; % launch eeglab

% Load data saved using Bittium Neurone EEG system with 64 channels of EEG and 3 channels of EMG data
EEG = pop_readneurone(fullfile('data_neurone', 'SPACETIME_020', 'NeurOne-2024-10-17T101023.ses'), 1);

% clock OK twice to accept standard channel locations for EEG sensors, EMG sensors will not have a location

BADTRIAL = 1:17; % first 17 trials contain artifacts from impedance measurement

% select EMG channels by name and create epochs
EMG = pop_select(EEG, 'channel', {'APBr', 'FDIr', 'ADMr'});
EMG = pop_epoch(EMG, {'A - Out'}, [-0.5 0.5]);
EMG = pop_select(EMG, 'notrial', BADTRIAL);
EMG = pop_rmbase(EMG, [-100 -10]);

% select EEG using channels having a location
EEG = pop_select(EEG, 'nochannel', find(cellfun(@isempty, {EEG.chanlocs.radius})));

EEG_epoched = pop_epoch(EEG, {'A - Out'}, [-1.5 1.5]);
EEG_epoched = pop_select(EEG_epoched, 'notrial', BADTRIAL);
EEG_epoched = pop_rmbase(EEG_epoched, [-1500 -10]);

% Save 5-channel EEGLAB EEG data for quality control

EEG = pop_select(EEG_epoched, 'channel', {'C3', 'FC1', 'FC5', 'CP1', 'CP5'});
save(fullfile('data_epoched', 'OPENLOOP_eeglab.mat'), 'EMG', 'EEG', '-v7.3')

% Create and save prestimulus fieldtrip data

EEG_prestim = pop_epoch(EEG_epoched, {'A - Out'}, [-1.5 -0.005]);
EEG_prestim = pop_resample(EEG_prestim, 1000); % this applies an anti-aliasing filter which will introduce small edge artifacts (the pre-stimulus data does not contain a problematic TMS artifact)

ftdata_emg = eeglab2fieldtrip(EMG, 'preprocessing');
ftdata_emg.trialinfo = [1:length(ftdata_emg.trial)]';
ftdata_emg.sampleinfo = [EMG.urevent(ismember({EMG.urevent.type}, {'A - Out'})).latency]' + ([-0.5 0.5]*EMG.srate);

ftdata_eeg_prestim = eeglab2fieldtrip(EEG_prestim, 'preprocessing');
ftdata_eeg_prestim.trialinfo = [1:length(ftdata_eeg_prestim.trial)]';
sampleinfo = [EEG_epoched.urevent(ismember({EEG_epoched.urevent.type}, {'A - Out'})).latency]' + ([-1.5 -0.005]*EEG_epoched.srate);
ftdata_eeg_prestim.sampleinfo = ceil(sampleinfo/5);

save(fullfile('data_epoched', 'OPENLOOP_ftdata.mat'), 'ftdata_emg', 'ftdata_eeg_prestim', '-v7.3')
