function eegplugin_fooof(fig,try_strings,catch_strings)
    % keyboard
    % tools menu
    eegfooofcommand = [try_strings.no_check '[EEG LASTCOM] = pop_eeg_fooof(EEG);' catch_strings.store_and_hist];
    
    tools_menu = findobj(fig, 'tag', 'tools');
    submenu = uimenu(tools_menu, 'Label', 'FOOOF', 'separator', 'on');
    uimenu(submenu, 'label', 'Fit power spectra', 'callback', eegfooofcommand); %then add eegh stuff
    uimenu(submenu, 'label', 'Plot fit', 'separator','on', 'callback', 'pop_eeg_fooofplot(EEG);');
    % study menu
    std_menu = findobj(fig, 'tag', 'study');
    submenu = uimenu(std_menu, 'Label', 'FOOOF', 'separator', 'on');
    uimenu( submenu, 'label', 'Fit power spectra', 'callback', 'pop_std_fooof(STUDY, ALLEEG)');
    uimenu(submenu, 'label', 'Plot cluster fit', 'separator','on', 'callback', 'pop_std_fooofplot(STUDY, ALLEEG);');