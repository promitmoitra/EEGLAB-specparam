% fooof_check_settings() - Check a struct of settings for the FOOOF model.
%
% Usage:
%  >> setting = fooof_check_settings(setting)
%
% Inputs:
%   setting        = struct, can optionally include:
%       setting.peak_width_limts
%       setting.max_n_peaks
%       setting.min_peak_height
%       setting.peak_threshold
%       setting.aperiodic_mode
%       setting.verbose
%
% Outputs:
%   setting        = struct, with all settings defined:
%       setting.peak_width_limts
%       setting.max_n_peaks
%       setting.min_peak_height
%       setting.peak_threshold
%       setting.aperiodic_mode
%       setting.verbose
%
% Notes:
%   This is a helper function, probably not called directly by the user.
%   Any settings not specified are set to default values

function setting = fooof_check_settings(setting)

    % Set defaults for all setting
    defaults = struct(...
        'peak_width_limits', [0.5, 12], ...
        'max_n_peaks', Inf, ...
        'min_peak_height', 0.0, ...
        'peak_threshold', 2.0, ...
        'aperiodic_mode', 'fixed', ...
        'verbose', true);

    % Overwrite any non-existent or nan setting with defaults
    for field = fieldnames(defaults)'
        if ~isfield(setting, field) || all(isnan(setting.(field{1})))
            setting.(field{1}) = defaults.(field{1});
        end
    end

end