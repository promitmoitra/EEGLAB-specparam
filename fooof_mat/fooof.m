% fooof() - Fit the FOOOF model on a neural power spectrum.
%
% Usage:
%   >> fooof_results = fooof(freqs, power_spectrum, f_range, setting, return_model);
%
% Inputs:
%   freqs           = row vector of frequency values
%   power_spectrum  = row vector of power values
%   f_range         = fitting range (Hz)
%   setting        = fooof model settings, in a struct, including:
%       setting.peak_width_limts
%       setting.max_n_peaks
%       setting.min_peak_height
%       setting.peak_threshold
%       setting.aperiodic_mode
%       setting.verbose
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
%     Not all settings need to be defined by the user.
%     Any settings that are not provided are set to default values.
%     To run with all defaults, input settings as an empty struct.

function fooof_results = fooof(freqs, power_spectrum, f_range, setting, return_model)
    % Check setting - get defaults for those not provided
    setting = fooof_check_settings(setting);

    % Convert inputs
    py_freqs = py.numpy.array(freqs);
    py_power_spectrum = py.numpy.array(power_spectrum);
    f_range = py.list(f_range);

    % Initialize FOOOF object
    fm = py.fooof.FOOOF(setting.peak_width_limits, ...
                        setting.max_n_peaks, ...
                        setting.min_peak_height, ...
                        setting.peak_threshold, ...
                        setting.aperiodic_mode, ...
                        setting.verbose);

    % Run FOOOF fit
    fm.fit(py_freqs, py_power_spectrum, f_range)

    % Extract outputs
    fooof_results = fm.get_results();
    fooof_results = fooof_unpack_results(fooof_results);
    offset = fooof_results.aperiodic_params(1); exponent = fooof_results.aperiodic_params(2);
    param_ap_fit = 10^offset*freqs.^(-exponent);fooof_results.ap_fit = param_ap_fit;
    py_ap_fit = py.numpy.array(param_ap_fit);

    % Adding band powers to fooof_results: # -> frequency range
    delta_band = py.list([1,4]);    %####
    theta_band = py.list([4,8]);    %####
    alpha_band = py.list([8,14]);   %######
    lbeta_band = py.list([14,20]);  %######
    hbeta_band = py.list([20,30]);  %##########
    lgamma_band = py.list([30,45]); %###############
    hgamma_band = py.list([45,80]); %###################################

    delta_psd = py.fooof.utils.trim_spectrum(py_freqs,py_power_spectrum,delta_band);
    theta_psd = py.fooof.utils.trim_spectrum(py_freqs,py_power_spectrum,theta_band);
    alpha_psd = py.fooof.utils.trim_spectrum(py_freqs,py_power_spectrum,alpha_band);
    lbeta_psd = py.fooof.utils.trim_spectrum(py_freqs,py_power_spectrum,lbeta_band);
    hbeta_psd = py.fooof.utils.trim_spectrum(py_freqs,py_power_spectrum,hbeta_band);
    lgamma_psd = py.fooof.utils.trim_spectrum(py_freqs,py_power_spectrum,lgamma_band);
    hgamma_psd = py.fooof.utils.trim_spectrum(py_freqs,py_power_spectrum,hgamma_band);
    
    fooof_results.bandpowers = ...
    [trapz(double(delta_psd{1}),double(delta_psd{2}));
     trapz(double(theta_psd{1}),double(theta_psd{2}));
     trapz(double(alpha_psd{1}),double(alpha_psd{2}));
     trapz(double(lbeta_psd{1}),double(lbeta_psd{2}));
     trapz(double(hbeta_psd{1}),double(hbeta_psd{2}));
     trapz(double(lgamma_psd{1}),double(lgamma_psd{2}));
     trapz(double(hgamma_psd{1}),double(hgamma_psd{2}))];

    delta_apsd = py.fooof.utils.trim_spectrum(py_freqs,py_ap_fit,delta_band);
    theta_apsd = py.fooof.utils.trim_spectrum(py_freqs,py_ap_fit,theta_band);
    alpha_apsd = py.fooof.utils.trim_spectrum(py_freqs,py_ap_fit,alpha_band);
    lbeta_apsd = py.fooof.utils.trim_spectrum(py_freqs,py_ap_fit,lbeta_band);
    hbeta_apsd = py.fooof.utils.trim_spectrum(py_freqs,py_ap_fit,hbeta_band);
    lgamma_apsd = py.fooof.utils.trim_spectrum(py_freqs,py_ap_fit,lgamma_band);
    hgamma_apsd = py.fooof.utils.trim_spectrum(py_freqs,py_ap_fit,hgamma_band);

    fooof_results.ap_bandpowers = ...
    [trapz(double(delta_apsd{1}),double(delta_apsd{2}));
     trapz(double(theta_apsd{1}),double(theta_apsd{2}));
     trapz(double(alpha_apsd{1}),double(alpha_apsd{2}));
     trapz(double(lbeta_apsd{1}),double(lbeta_apsd{2}));
     trapz(double(hbeta_apsd{1}),double(hbeta_apsd{2}));
     trapz(double(lgamma_apsd{1}),double(lgamma_apsd{2}));
     trapz(double(hgamma_apsd{1}),double(hgamma_apsd{2}))];

%     fooof_results.periodic_bandpowers = ...
%     [trapz(double(delta_psd{1}),double(delta_psd{2}))-trapz(double(delta_apsd{1}),double(delta_apsd{2}));
%      trapz(double(theta_psd{1}),double(theta_psd{2}))-trapz(double(theta_apsd{1}),double(theta_apsd{2}));
%      trapz(double(alpha_psd{1}),double(alpha_psd{2}))-trapz(double(alpha_apsd{1}),double(alpha_apsd{2}));
%      trapz(double(lbeta_psd{1}),double(lbeta_psd{2}))-trapz(double(lbeta_apsd{1}),double(lbeta_apsd{2}));
%      trapz(double(hbeta_psd{1}),double(hbeta_psd{2}))-trapz(double(hbeta_apsd{1}),double(hbeta_apsd{2}));
%      trapz(double(lgamma_psd{1}),double(lgamma_psd{2}))-trapz(double(lgamma_apsd{1}),double(lgamma_apsd{2}))];

    % Re-calculating r-squared
    %   r_squared doesn't seem to get computed properly (in NaN).
    %   It is unclear why this happens, other than the error can be traced
    %   back to the internal call to `np.cov`, and fails when this function
    %   gets two arrays as input.
    %   Therefore, we can simply recalculate r-squared
    coefs = corrcoef(double(py.array.array('d', fm.power_spectrum)), ...
                     double(py.array.array('d', fm.fooofed_spectrum_)));
    fooof_results.r_squared = coefs(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Also return the actual model fit, if requested
    % This will default to not return model, if variable not set

%%% Generates incorrect spectra with negative abs power density (uV^2/Hz) values :( 
%%% even though the estimated aperiodic_parameters are correct! Commenting for now...    
%     if exist('return_model', 'var') && return_model
% %         % Get the model, and add outputs to fooof_results
% %         model_out = fooof_get_model(fm);
% %         for field = fieldnames(model_out)'
% %             fooof_results.(field{1}) = model_out.(field{1});
% %         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
