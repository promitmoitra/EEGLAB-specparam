% fooof() - Fit the FOOOF model on a neural power spectrum.
%
% Usage:
%   >> fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
%
% Inputs:
%   freqs           = row vector of frequency values
%   power_spectrum  = row vector of power values
%   f_range         = fitting range (Hz)
%   settings        = fooof model settings, in a struct, including:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_height
%       settings.peak_threshold
%       settings.aperiodic_mode
%       settings.verbose
%   return_model    = boolean of whether to return the FOOOF model fit, optional
%
% Outputs:
%   fooof_results   = fooof model ouputs, in a struct, including:
%       fooof_results.aperiodic_params
%       fooof_results.peak_params
%       fooof_results.gaussian_params
%       fooof_results.error
%       fooof_results.r_squared
%       if return_model is true, it also includes:
%            fooof_results.freqs
%            fooof_results.power_spectrum
%            fooof_results.fooofed_spectrum
%            fooof_results.ap_fit
%
% Notes
%   Not all settings need to be defined by the user.
%     Any settings that are not provided are set to default values.
%     To run with all defaults, input settings as an empty struct.

function fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model)

    % Check settings - get defaults for those not provided
    settings = fooof_check_settings(settings);

    % Convert inputs
    freqs = py.numpy.array(freqs);
    power_spectrum = py.numpy.array(power_spectrum);
    f_range = py.list(f_range);

    % Initialize FOOOF object
    fm = py.fooof.FOOOF(settings.peak_width_limits, ...
                        settings.max_n_peaks, ...
                        settings.min_peak_height, ...
                        settings.peak_threshold, ...
                        settings.aperiodic_mode, ...
                        settings.verbose);

    % Run FOOOF fit
    fm.fit(freqs, power_spectrum, f_range)

    % Extract outputs
    fooof_results = fm.get_results();
    fooof_results = fooof_unpack_results(fooof_results);
    
    % Adding band powers to fooof_results:
    delta_band = py.list([1,4]);delta_power = py.fooof.utils.trim_spectrum(freqs,power_spectrum,delta_band);
    theta_band = py.list([4,7]);theta_power = py.fooof.utils.trim_spectrum(freqs,power_spectrum,theta_band);
    alpha_band = py.list([7,14]);alpha_power = py.fooof.utils.trim_spectrum(freqs,power_spectrum,alpha_band);
    beta_band = py.list([15,30]);beta_power = py.fooof.utils.trim_spectrum(freqs,power_spectrum,beta_band);
    gamma_band = py.list([30,45]);gamma_power = py.fooof.utils.trim_spectrum(freqs,power_spectrum,gamma_band);
%     fooof_results.bandfreqs =  {double(delta_power{1});
%                                 double(theta_power{1});
%                                 double(alpha_power{1});
%                                 double(beta_power{1});
%                                 double(gamma_power{1})};
    fooof_results.bandpowers = [mean(double(delta_power{2}));
                                mean(double(theta_power{2}));
                                mean(double(alpha_power{2}));
                                mean(double(beta_power{2}));
                                mean(double(gamma_power{2}))];

    % Re-calculating r-squared
    %   r_squared doesn't seem to get computed properly (in NaN).
    %   It is unclear why this happens, other than the error can be traced
    %   back to the internal call to `np.cov`, and fails when this function
    %   gets two arrays as input.
    %   Therefore, we can simply recalculate r-squared
    coefs = corrcoef(double(py.array.array('d', fm.power_spectrum)), ...
                     double(py.array.array('d', fm.fooofed_spectrum_)));
    fooof_results.r_squared = coefs(2);

    % Also return the actual model fit, if requested
    %   This will default to not return model, if variable not set
    if exist('return_model', 'var') && return_model

        % Get the model, and add outputs to fooof_results
        model_out = fooof_get_model(fm);
        for field = fieldnames(model_out)'
            fooof_results.(field{1}) = model_out.(field{1});
        end

    end

end
