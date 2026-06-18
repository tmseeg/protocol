%% Extract MEP amplitude

% Note that this step is not part of the protocol but is included here for illustrative purposes

assert(size(ftdata_emg.trialinfo,1) == size(ftdata_eeg_prestim.trialinfo,1), 'eeg and emg trialinfo indices do not match');

datacube_emg = cat(3, ftdata_emg.trial{:});

MEP_RANGE = 2620:2680; % these are the trial indices for determining maximum and minimum MEP amplitude in this participant
CHANNEL = 2;

mep_amplitude = squeeze(range(datacube_emg(CHANNEL, MEP_RANGE, :), 2));

clear('datacube_emg', 'MEP_RANGE', 'CHANNEL')

save('mep_amplitude.mat', 'mep_amplitude')