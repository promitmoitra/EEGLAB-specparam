% fooof_group() - Run the fooof model on a group of neural power spectra.
%
% Usage:
%   fooof_results = fooof_group(freqs, psds, f_range, setting);
%
% Inputs:
%   freqs           = row vector of frequency values
%   psds            = matrix of power values, which each row representing a spectrum
%   f_range         = fitting range (Hz)
%   setting        = fooof model setting, in a struct, including:
%       setting.peak_width_limts
%       setting.max_n_peaks
%       setting.min_peak_height
%       setting.peak_threshold
%       setting.aperiodic_mode
%       setting.verbose
%
% Outputs:
%   fooof_results   = fooof model ouputs, in a struct, including:
%       fooof_results.aperiodic_params
%       fooof_results.peak_params
%       fooof_results.gaussian_params
%       fooof_results.error
%       fooof_results.r_squared
%
% Notes
%   Not all setting need to be set. Any setting that are not
%     provided as set to default values. To run with all defaults,
%     input setting as an empty struct.

function fooof_results = fooof_group(freqs, psds, f_range, setting, return_model)
    if ~exist('return_model', 'var')
        return_model = false; 
    end

    % Check setting - get defaults for those not provided
    setting = fooof_check_settings(setting);

    % Initialize object to collect FOOOF results
    fooof_results = [];

    % Run FOOOF across the group of power spectra
    for psd = psds
%         sprintf('fooof_group: freqs %d,',size(freqs))
%         sprintf('fooof_group: psd` %d,',size(psd'))
        cur_results = fooof(freqs, psd', f_range, setting, return_model);
        fooof_results = [fooof_results, cur_results];
    end

end