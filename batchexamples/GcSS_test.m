
% select localizer SPM.mat files
subjects = {...
    '/mindhive/evlab/u/Shared/SUBJECTS/832_FED_20200220c_3T2_PL2017/firstlevel_langlocSN/SPM.mat',...
    '/mindhive/evlab/u/Shared/SUBJECTS/833_FED_20200224a_3T2_PL2017/firstlevel_langlocSN/SPM.mat',...
    '/mindhive/evlab/u/Shared/SUBJECTS/834_FED_20200303a_3T2_PL2017/firstlevel_langlocSN/SPM.mat',...
    '/mindhive/evlab/u/Shared/SUBJECTS/835_FED_20200305a_3T2_PL2017/firstlevel_langlocSN/SPM.mat',...
    '/mindhive/evlab/u/Shared/SUBJECTS/836_FED_20210308a_3T2_PL2017/firstlevel_langlocSN/SPM.mat',...
    '/mindhive/evlab/u/Shared/SUBJECTS/837_FED_20210629a_3T1_PL2017/firstlevel_langlocSN/SPM.mat',...
};

% setup SPM_SS control structure
ss=struct(...
    'swd',fullfile(pwd,'GcSS_test'),... % working directory
    'files_spm',{subjects},... % SPM structures for subjects
    'EffectOfInterest_contrasts',{{'S','N'}},... % name of effect contrasts in SPM struct
    'Localizer_contrasts',{{'S-N'}},... % name of localizer contrast in SPM struct
    'Localizer_thr_type',{{'none'}},... % how to select voxels, e.g. sig threshold, percentile, etc., in this case, uncorrected sig threshold
    'Localizer_thr_p',0.001,...  % cutoff defined here, can be p-val, percentile, nvoxels, etc. based on specification above
    'overlap_thr_roi',0,... % minimum proportion of subjects showing a significant localizer effect within an ROI when measuring effect sizes; 0 means don't disregard any ROIs.
    'overlap_thr_vox',0.1,... % minimum proportion of subjects showing an effect in a voxel when constructing ROI parcellation.
    'type','GcSS',... % analysis type
    'smooth',8,... % smoothing in Gaussian FWHM mm
    'overwrite', 1,... % overwrite output files if run again
    'model',1,... % one-sample t-test (for single-group analyses)
    'ExplicitMasking', [], ... % if you want to apply additional explicit masking to constrain GcSS
    'estimation','OLS'... % ordinary least squares estimation
);

% run GcSS
addpath('/om/group/evlab/software/spm12');
addpath('/om/group/evlab/software/spm_ss');
ss=spm_ss_design(ss);
ss=spm_ss_estimate(ss);
