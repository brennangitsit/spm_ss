function ss=spm_ss_estimate_ROI(ss)
% SPM_SS_ESTIMATE_ROI subject-specific ROI-based model estimation
%
% ss=spm_ss_estimate_ROI(ss)
% see SPM_SS_DESIGN, SPM_SS_ESTIMATE

%% Load analysis file if not provided
if nargin<1,
    str='Select spm_ss*.mat analysis file';
    disp(str);
    Pdefault='';objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_ss'),Pdefault=objdata.files_spm_ss;end;end;
    P=spm_select(1,'^SPM_ss.*\.mat$',str,{Pdefault});
    if numel(objname)==1&&~isempty(P),objdata.files_spm_ss=P;set(objname,'userdata',objdata);end;
    load(P);
    ss.swd=fileparts(P);
end

cwd=pwd;

%% Load explicit mask (optional additional constraint applied to all subjects)
if ~isempty(ss.ExplicitMasking),XM=spm_vol(ss.ExplicitMasking);else XM=[];end

%% Build file lists for effect-of-interest and localizer volumes
% ss.PY: cell array of effect-of-interest contrast file paths
%        ordered as [effects x subjects x cross-validation folds]
% ss.PN: cell array of localizer mask file paths
% ss.PV: [nsubjects x total_folds] weight matrix (1/nfolds per subject)
neffects=size(ss.EffectOfInterest{1},1);
k=0;
for ne=1:neffects,
    for n=1:ss.n,
        for m=1:numel(ss.Localizer{n}),
            [pth2,nm2,ext2,num2]=spm_fileparts(ss.EffectOfInterest{n}{ne,m});
            Yvolume=[nm2,ext2,num2];
            k=k+1;
            ss.PY{k}=fullfile(pth2,Yvolume);
        end
    end
end
k=0;
ss.PV=[];
ss.PV_subjnames={};
for n=1:ss.n,
    for m=1:numel(ss.Localizer{n}),
        [pth1,nm1,ext1,num1]=spm_fileparts(ss.Localizer{n}{m});
        Nvolume=[nm1,ext1,num1];
        k=k+1;
        ss.PN{k}=fullfile(pth1,Nvolume);
        ss.PV(n,k)=1/numel(ss.Localizer{n});  % equal weight per fold
    end
    ss.PV_subjnames{n}=num2str(n);
    if numel(ss.Localizer{n})>0
        [pth1a,pth1b]=fileparts(pth1);
        if ~isempty(regexp(pth1b,'^model_|^firstlevel_'))&&~isempty(pth1a), [pth1b,pth1a]=fileparts(pth1a); ss.PV_subjnames{n}=pth1a; end
    end
end

%% Load volume headers
cd(ss.swd);
ss.VN=spm_vol(char(ss.PN));                             % localizer masks: [total_folds_across_subjects x 1]
ss.VY=reshape(spm_vol(char(ss.PY)),numel(ss.VN),[]);    % effects: [total_folds_across_subjects x neffects]
ss.PY=reshape(ss.PY,numel(ss.VN),[]);

%% Check spatial alignment across subjects
% refspace tracks whether effect volumes share the same space as each other
% and as the localizer mask, to optimise data reading later
ss.refspace.masks={};
ss.refspace.effects={};
ss.refspace.effectsinsamespace=false(1,ss.n);
ss.refspace.effectsinsamespaceasmask=false(1,ss.n);
for n=1:ss.n,
    idxk=ss.PV(n,:)>0;
    idxk1=find(idxk,1);
    if ~isempty(idxk1)
        [pth1,nm1,ext1,num1]=spm_fileparts(ss.VY(idxk1,1).fname);
        [pth1a,pth1b]=fileparts(pth1);
        if ~isempty(regexp(pth1b,'^model_|^firstlevel_'))&&~isempty(pth1a), [pth1b,pth1a]=fileparts(pth1a); ss.PV_subjnames{n}=pth1a; end
        nvy=ss.VN(idxk1);
        ss.refspace.masks{n}=struct('mat',nvy.mat,'dim',nvy.dim);
        nvy=ss.VY(idxk,:);
        ss.refspace.effects{n}=struct('mat',{nvy.mat},'dim',{nvy.dim});
        try, ss.refspace.effectsinsamespace(n)=all(reshape(~diff(cat(4,nvy.mat),1,4),[],1))&all(reshape(~diff(cat(4,nvy.dim),1,4),[],1)); end
        try, ss.refspace.effectsinsamespaceasmask(n)=all(reshape(~diff(cat(4,ss.VN(idxk1).mat,nvy.mat),1,4),[],1))&all(reshape(~diff(cat(4,ss.VN(idxk1).dim,nvy.dim),1,4),[],1)); end
    end
end

%% Determine cross-validation structure
% ss.crossvalidation: 'kfold', 'loso', or 'firstsess'
if ~isfield(ss, 'crossvalidation'), ss.crossvalidation = 'kfold'; end

% Determine which mode to use
use_loso = (ss.typen == 2) && strcmpi(ss.crossvalidation, 'loso') && (ss.n > 2);
use_firstsess = (ss.typen == 2) && strcmpi(ss.crossvalidation, 'firstsess');

if use_firstsess
    fprintf('Cross-validation: firstsess (parcellation/localization uses sessions 2+, effects from session 1 only)\n');
elseif use_loso
    fprintf('Cross-validation: loso (each subject gets parcels computed from N-1 other subjects)\n');
else
    fprintf('Cross-validation: kfold (all subjects use the same parcels; k-fold cross-validation for fROIs)\n');
end

%% Create inter-subject overlap maps
% For standard mode: single overlap map using all subjects
% For LOSO mode: N+1 overlap maps (one leaving out each subject, plus the full map)
% For first-session holdout: use only fold 1's localizer (ORTH_TO_SESSION01) per subject
ssPM1 = ['Overlap', ext1];

% Read subject localizers into memory
% - Standard/LOSO: average across all folds per subject
% - First-session holdout: use only fold 1 (ORTH_TO_SESSION01, i.e., sessions 2+)
subj_localizers = cell(1, ss.n);
for n = 1:ss.n
    idx = find(ss.PV(n,:));

    if use_firstsess
        % First-session holdout: only use fold 1's localizer (ORTH_TO_SESSION01)
        k = 1;  % only read the first fold
        a1 = ss.VN(idx(k));
        if n == 1
            b1 = spm_read_vols(a1);
            ref_mat = a1.mat;
            ref_dim = a1.dim;
            [XYZ{1}, XYZ{2}, XYZ{3}] = ndgrid(1:a1.dim(1), 1:a1.dim(2), 1:a1.dim(3));
            tXYZ = reshape(cat(4, XYZ{:}), [], 3)';
            tXYZ = a1.mat * cat(1, tXYZ, ones(1, size(tXYZ, 2)));
        else
            b1 = reshape(spm_get_data(a1, pinv(a1.mat) * tXYZ), ref_dim);
        end
        p1 = (b1 > 0);  % binary localizer from fold 1
    else
        % Standard/LOSO: average across all folds
        p1 = 0;
        for k = 1:numel(idx)
            a1 = ss.VN(idx(k));
            if n == 1 && k == 1
                b1 = spm_read_vols(a1);
                ref_mat = a1.mat;
                ref_dim = a1.dim;
                [XYZ{1}, XYZ{2}, XYZ{3}] = ndgrid(1:a1.dim(1), 1:a1.dim(2), 1:a1.dim(3));
                tXYZ = reshape(cat(4, XYZ{:}), [], 3)';
                tXYZ = a1.mat * cat(1, tXYZ, ones(1, size(tXYZ, 2)));
            elseif n == 1
                b1 = spm_read_vols(a1);
            else
                b1 = reshape(spm_get_data(a1, pinv(a1.mat) * tXYZ), ref_dim);
            end
            p1 = p1 + (b1 > 0) / numel(idx);
        end
    end
    subj_localizers{n} = p1;  % localizer for subject n
end

% Compute the full overlap map (all subjects)
p_full = zeros(ref_dim);
for n = 1:ss.n
    p_full = p_full + (subj_localizers{n} > 0.5);  % binarize at 0.5
end
overlap_full = p_full / ss.n;

% Save the full overlap map
e1 = struct('fname', ssPM1, 'descrip', 'spm_ss (inter-subject overlap map)', ...
    'mat', ref_mat, 'dim', ref_dim, 'dt', [spm_type('float32'), spm_platform('bigend')]);
e1 = spm_write_vol(e1, overlap_full);

% Create per-subject localizer QA volume
ss.PL = ['Localizer', ext1];
e0 = struct('fname', ss.PL, 'descrip', 'spm_ss (localizer mask for each subject)', ...
    'mat', ref_mat, 'dim', ref_dim, 'n', [1,1], 'pinfo', [1;0;0], ...
    'dt', [spm_type('float32'), spm_platform('bigend')]);
try, spm_unlink(e0.fname); end
e0 = repmat(e0, [ss.n, 1]);
for nb = 1:ss.n, e0(nb).n = [nb, 1]; end
e0 = spm_create_vol(e0);
for n = 1:ss.n
    spm_write_vol(e0(n), subj_localizers{n});
end

% For LOSO: compute leave-one-out overlap maps
if use_loso
    overlap_loso = cell(1, ss.n);
    for left_out = 1:ss.n
        p_loso = zeros(ref_dim);
        for n = 1:ss.n
            if n ~= left_out
                p_loso = p_loso + (subj_localizers{n} > 0.5);
            end
        end
        overlap_loso{left_out} = p_loso / (ss.n - 1);
    end
    % Save LOSO overlap maps for QA
    for left_out = 1:ss.n
        e_loso = struct('fname', sprintf('Overlap_loso%02d%s', left_out, ext1), ...
            'descrip', sprintf('spm_ss (overlap map leaving out subject %d)', left_out), ...
            'mat', ref_mat, 'dim', ref_dim, 'dt', [spm_type('float32'), spm_platform('bigend')]);
        spm_write_vol(e_loso, overlap_loso{left_out});
    end
end

%% Homologous parcellation: make overlap maps bilaterally symmetric
% When ss.homologous is set, flip the overlap map(s) across the sagittal
% plane (dimension 1 = left-right) and average with the original.
if ss.typen == 2 && isfield(ss, 'homologous') && ss.homologous
    disp('Creating homologous (bilaterally symmetric) overlap map(s)...');
    % Transform full overlap
    overlap_data = overlap_full;
    flipped_data = flip(overlap_data, 1);
    overlap_full = (overlap_data + flipped_data) / 2;
    spm_write_vol(e1, overlap_full);

    % Transform LOSO overlaps if applicable
    if use_loso
        for left_out = 1:ss.n
            overlap_data = overlap_loso{left_out};
            flipped_data = flip(overlap_data, 1);
            overlap_loso{left_out} = (overlap_data + flipped_data) / 2;
        end
    end
end

%% Smooth and watershed to create parcels (GcSS only)
ss.PM1 = ssPM1;
parcels_loso = cell(1, ss.n);  % LOSO parcels (one per subject)

if ss.typen == 2
    % Prepare explicit mask in reference space
    if ~isempty(XM)
        [XYZ{1}, XYZ{2}, XYZ{3}] = ndgrid(1:ref_dim(1), 1:ref_dim(2), 1:ref_dim(3));
        tXYZ_mask = reshape(cat(4, XYZ{:}), [], 3)';
        tXYZ_mask = cat(1, tXYZ_mask, ones(1, size(tXYZ_mask, 2)));
        tXM = reshape(double(spm_get_data(XM, pinv(XM.mat) * ref_mat * tXYZ_mask) > 0), ref_dim(1:3));
    else
        tXM = 1;
    end

    % Create the main (full) parcellation
    ss.PM2 = ['sOverlap', ext1];
    spm_smooth(ssPM1, ss.PM2, ss.smooth * [1, 1, 1]);
    a2 = spm_vol(ss.PM2);
    b2 = spm_read_vols(a2);
    b3 = spm_ss_watershed(-b2, find(b2 >= ss.overlap_thr_vox & tXM), ss.boundary_roi, ss.minsize_roi);
    ss.PM = ['fROIs', ext1];
    a3 = struct('fname', ss.PM, 'mat', a2.mat, 'dim', a2.dim, ...
        'dt', [spm_type('int16'), spm_platform('bigend')], 'pinfo', [1;0;0]);
    spm_write_vol(a3, b3);

    if isfield(ss, 'homologous') && ss.homologous
        fprintf('GcSS defined %d homologous ROIs (full parcellation)\n', max(b3(:)));
    else
        fprintf('GcSS defined %d ROIs (full parcellation)\n', max(b3(:)));
    end

    % Create LOSO parcellations if enabled
    if use_loso
        fprintf('Creating LOSO parcellations (one per subject)...\n');
        for left_out = 1:ss.n
            % Write LOSO overlap
            loso_overlap_fname = sprintf('Overlap_loso%02d%s', left_out, ext1);
            e_loso = struct('fname', loso_overlap_fname, 'mat', ref_mat, 'dim', ref_dim, ...
                'dt', [spm_type('float32'), spm_platform('bigend')]);
            spm_write_vol(e_loso, overlap_loso{left_out});

            % Smooth
            loso_soverlap_fname = sprintf('sOverlap_loso%02d%s', left_out, ext1);
            spm_smooth(loso_overlap_fname, loso_soverlap_fname, ss.smooth * [1, 1, 1]);

            % Watershed
            a2_loso = spm_vol(loso_soverlap_fname);
            b2_loso = spm_read_vols(a2_loso);
            b3_loso = spm_ss_watershed(-b2_loso, find(b2_loso >= ss.overlap_thr_vox & tXM), ss.boundary_roi, ss.minsize_roi);
            parcels_loso{left_out} = b3_loso;

            % Save LOSO parcels
            loso_parcel_fname = sprintf('fROIs_loso%02d%s', left_out, ext1);
            a3_loso = struct('fname', loso_parcel_fname, 'mat', a2_loso.mat, 'dim', a2_loso.dim, ...
                'dt', [spm_type('int16'), spm_platform('bigend')], 'pinfo', [1;0;0]);
            spm_write_vol(a3_loso, b3_loso);

            fprintf('  Subject %d LOSO: %d ROIs\n', left_out, max(b3_loso(:)));
        end
    end
else
    % Manual ROI mode: use user-provided ROI file(s)
    ss.PM = ss.ManualROIs;
    use_loso = false;
end
ss.VM = spm_vol(char(ss.PM));
fprintf(1, '\n');

%% Map fROIs into each subject's native space and create QA images
% For LOSO parcellation:
%   frois{n} = full parcellation mapped to subject n's space (for QA)
%   frois_loso{n} = LOSO parcellation for subject n (computed without subject n)
% For standard parcellation:
%   frois{n} = same parcellation for all subjects

if use_loso
    frois_loso = cell(1, ss.n);
end

for n = 1:ss.n
    % resample fROI labels into subject n's localizer mask space
    idxk = ss.PV(n,:) > 0;
    idxk1 = find(idxk, 1);
    [XYZ{1}, XYZ{2}, XYZ{3}] = ndgrid(1:ss.VN(idxk1).dim(1), 1:ss.VN(idxk1).dim(2), 1:ss.VN(idxk1).dim(3));
    tXYZ = reshape(cat(4, XYZ{:}), [], 3)';
    tXYZ = cat(1, tXYZ, ones(1, size(tXYZ, 2)));

    % Map full parcellation (for QA and backward compatibility)
    vm = ss.VM(min(numel(ss.VM), n));
    frois{n} = reshape(round(spm_get_data(vm, pinv(vm.mat) * ss.VN(idxk1).mat * tXYZ)), ss.VN(idxk1).dim(1:3));
    if n == 1, nrois = max(frois{n}(:));
    else nrois = max(nrois, max(frois{n}(:)));
    end

    % Map LOSO parcellation for this subject (if enabled)
    if use_loso
        loso_parcel_fname = sprintf('fROIs_loso%02d%s', n, ext1);
        vm_loso = spm_vol(loso_parcel_fname);
        frois_loso{n} = reshape(round(spm_get_data(vm_loso, pinv(vm_loso.mat) * ss.VN(idxk1).mat * tXYZ)), ss.VN(idxk1).dim(1:3));
    end

    % zero out ROI labels outside explicit mask
    if ~isempty(XM)
        tXM_subj = reshape(double(spm_get_data(XM, pinv(XM.mat) * ss.VN(idxk1).mat * tXYZ) > 0), ss.VN(idxk1).dim(1:3));
        frois{n} = tXM_subj .* frois{n};
        if use_loso
            frois_loso{n} = tXM_subj .* frois_loso{n};
        end
    end

    % QA image: parcel labels in subject space (using full parcellation for visualization)
    a4 = struct('fname', [sprintf('QA_parcels.%s', ss.PV_subjnames{n}), ext1], ...
        'mat', ss.refspace.masks{n}.mat, 'dim', ss.refspace.masks{n}.dim, ...
        'dt', [spm_type('float32'), spm_platform('bigend')], ...
        'descrip', sprintf('QA display parcels for subject#%d', n), 'pinfo', [1;0;0]);
    spm_write_vol(a4, reshape(frois{n}, ss.VN(idxk1).dim(1:3)));

    % QA image: average effect across all contrasts
    a4 = struct('fname', [sprintf('QA_effects.%s', ss.PV_subjnames{n}), ext1], ...
        'mat', ss.refspace.masks{n}.mat, 'dim', ss.refspace.masks{n}.dim, ...
        'dt', [spm_type('float32'), spm_platform('bigend')], ...
        'descrip', sprintf('QA display average all effects for subject#%d', n), 'pinfo', [1;0;0]);
    tvy = reshape(ss.VY(idxk,:), [], 1);
    if isempty(idxk1)
        tY = zeros(numel(tvy), size(tXYZ, 2));
    else
        if ss.refspace.effectsinsamespaceasmask(n)
            tY = spm_get_data(tvy, tXYZ);
        elseif ss.refspace.effectsinsamespace(n)
            tY = spm_get_data(tvy, pinv(tvy(1).mat) * ss.refspace.masks{n}.mat * tXYZ);
        else
            tY = zeros(numel(tvy), size(tXYZ, 2));
            for nk = 1:numel(tvy)
                tY(nk,:) = spm_get_data(tvy(nk), pinv(tvy(nk).mat) * ss.refspace.masks{n}.mat * tXYZ);
            end
        end
    end
    tY(isnan(tY)) = 0;
    spm_write_vol(a4, reshape(mean(tY, 1), ss.VN(idxk1).dim(1:3)));

    % QA image: average localizer mask
    a4 = struct('fname', [sprintf('QA_localizer.%s', ss.PV_subjnames{n}), ext1], ...
        'mat', ss.refspace.masks{n}.mat, 'dim', ss.refspace.masks{n}.dim, ...
        'dt', [spm_type('float32'), spm_platform('bigend')], ...
        'descrip', sprintf('QA display average localizer mask for subject#%d', n), 'pinfo', [1;0;0]);
    tvy = reshape(ss.VN(idxk), [], 1);
    if isempty(idxk1)
        tN = zeros(numel(tvy), size(tXYZ, 2));
    else
        tN = spm_get_data(tvy, tXYZ);
    end
    tN(isnan(tN)) = 0;
    tN = reshape(mean(tN, 1), ss.VN(idxk1).dim(1:3));
    spm_write_vol(a4, tN);
end

% Load ROI names/IDs from label file if available, otherwise use numeric labels
[ss.VM_roinames, ss.VM_roiids] = spm_ss_roilabels(ss.VM(1).fname);
if isempty(ss.VM_roinames)
    ss.VM_roinames = arrayfun(@num2str, 1:nrois, 'uni', 0);
    ss.VM_roiids = 1:nrois;
else
    nrois = numel(ss.VM_roiids);
end

%% GLM estimation: extract effect sizes per ROI
% Nb(1) = number of between-subjects regressors, Nb(2) = number of effects
Nb=[size(ss.X,2),neffects];
extname=['_',ss.type];

% initialise output volumes for beta estimates and overlap
VB=struct('fname',['spm_ss',extname,'_beta.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'n',[1,1],...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (effect sizes parameter estimates)');
try, spm_unlink(VB.fname); end
VB=repmat(VB,[prod(Nb),1]);for nb=1:prod(Nb),VB(nb).n=[nb,1];end
VB=spm_create_vol(VB);
VO=struct('fname',['spm_ss',extname,'_overlap.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (proportion overlap)');
VO=spm_create_vol(VO);

% pre-allocate per-ROI result matrices
% Bplane: beta estimates [nregressors x neffects x nrois]
% Cplane: whitening weights [nsubjects x nrois]
% Eplane: residual sum-of-squares [neffects x neffects x nrois]
% Oplane: inter-subject overlap proportion [1 x nrois]
% Zplane: per-subject effect sizes [nsubjects x neffects x nrois]
% Nplane: ROI voxel counts [nsubjects x nrois]
% Pplane: mean localizer coverage [1 x nrois]
% Mplane: per-fold localizer mask sizes [total_folds x nrois]
Bplane=nan+zeros([Nb,nrois]);Cplane=zeros(ss.n,nrois);Eplane=nan+zeros(Nb(2),Nb(2),nrois);Oplane=nan+zeros(1,nrois);Zplane=zeros(ss.n,Nb(2),nrois);Nplane=nan+zeros([ss.n,nrois]);Pplane=nan+zeros([1,nrois]);Mplane=nan+zeros([numel(ss.VN),nrois]);

fprintf('Performing model estimation');
if use_firstsess
    fprintf(' (first-session holdout mode)...');
elseif use_loso
    fprintf(' (with LOSO cross-validated parcellation)...');
else
    fprintf('...');
end

for nroi = 1:nrois
    fprintf(1, '.');

    % For first-session holdout: only use fold 1 (SESSION01 effects, ORTH_TO_SESSION01 localizer)
    % For standard/LOSO: use all folds and average
    if use_firstsess
        Y_subj = zeros(ss.n, Nb(2));
        N_subj = zeros(ss.n, 1);

        for n = 1:ss.n
            idxk = ss.PV(n,:) > 0;
            fold_indices = find(idxk);
            global_fold_idx = fold_indices(1);  % only use fold 1

            % Choose parcellation (LOSO disabled in first-session holdout mode)
            current_froi = frois{n};

            % Find voxels belonging to this ROI
            idx = find(current_froi == ss.VM_roiids(nroi));
            [idx1, idx2, idx3] = ind2sub(size(current_froi), idx);
            xyz = [idx1, idx2, idx3, ones(numel(idx1), 1)]';
            Nplane(n, nroi) = numel(idx);

            % Read localizer weight for fold 1 at ROI voxels
            tN_single = spm_get_data(ss.VN(global_fold_idx), xyz);

            % Read effect data for fold 1 (SESSION01) at ROI voxels
            tvy = ss.VY(global_fold_idx, :);
            if isempty(idx)
                tY_single = zeros(Nb(2), 0);
            else
                if ss.refspace.effectsinsamespaceasmask(n)
                    tY_single = spm_get_data(tvy(:), xyz);
                elseif ss.refspace.effectsinsamespace(n)
                    tY_single = spm_get_data(tvy(:), pinv(tvy(1).mat) * ss.refspace.masks{n}.mat * xyz);
                else
                    tY_single = zeros(Nb(2), size(xyz, 2));
                    for ne = 1:Nb(2)
                        tY_single(ne, :) = spm_get_data(tvy(ne), pinv(tvy(ne).mat) * ss.refspace.masks{n}.mat * xyz);
                    end
                end
            end

            % Zero out voxels where localizer mask is empty or NaN
            Z_mask = (tN_single == 0) | isnan(tN_single);
            tY_single(:, Z_mask) = 0;
            tN_single(Z_mask) = 0;
            tY_single(isnan(tY_single)) = 0;

            Mplane(global_fold_idx, nroi) = sum(tN_single, 2);

            % Weighted average across voxels within ROI
            if ~isempty(idx) && any(tN_single > 0)
                tY_weighted = mean(tY_single .* repmat(tN_single, [Nb(2), 1]), 2);
                tN_mean = mean(tN_single, 2);
                % Normalise by localizer weights
                Y_subj(n, :) = (tY_weighted / max(eps, tN_mean))';
                N_subj(n) = tN_mean;
            else
                Y_subj(n, :) = 0;
                N_subj(n) = 0;
            end
        end

        Y = Y_subj;
        N = N_subj;
        sN = mean(N > 1e-4, 1);

    else
        % Standard k-fold or LOSO mode: use all folds and average
        Y = zeros(numel(ss.VY), 1);
        N = zeros(numel(ss.VN), 1);

        for n = 1:ss.n
            idxk = ss.PV(n,:) > 0;
            idxk1 = find(idxk, 1);
            fold_indices = find(idxk);

            % Choose parcellation
            if use_loso
                current_froi = frois_loso{n};
            else
                current_froi = frois{n};
            end

            % Find voxels belonging to this ROI
            idx = find(current_froi == ss.VM_roiids(nroi));
            [idx1, idx2, idx3] = ind2sub(size(current_froi), idx);
            xyz = [idx1, idx2, idx3, ones(numel(idx1), 1)]';
            Nplane(n, nroi) = numel(idx);

            for local_fold = 1:numel(fold_indices)
                global_fold_idx = fold_indices(local_fold);

                tN_single = spm_get_data(ss.VN(global_fold_idx), xyz);
                tvy = ss.VY(global_fold_idx, :);
                if isempty(idx)
                    tY_single = zeros(Nb(2), 0);
                else
                    if ss.refspace.effectsinsamespaceasmask(n)
                        tY_single = spm_get_data(tvy(:), xyz);
                    elseif ss.refspace.effectsinsamespace(n)
                        tY_single = spm_get_data(tvy(:), pinv(tvy(1).mat) * ss.refspace.masks{n}.mat * xyz);
                    else
                        tY_single = zeros(Nb(2), size(xyz, 2));
                        for ne = 1:Nb(2)
                            tY_single(ne, :) = spm_get_data(tvy(ne), pinv(tvy(ne).mat) * ss.refspace.masks{n}.mat * xyz);
                        end
                    end
                end

                Z_mask = (tN_single == 0) | isnan(tN_single);
                tY_single(:, Z_mask) = 0;
                tN_single(Z_mask) = 0;
                tY_single(isnan(tY_single)) = 0;

                Mplane(global_fold_idx, nroi) = sum(tN_single, 2);

                if ~isempty(idx) && any(tN_single > 0)
                    tY_weighted = mean(tY_single .* repmat(tN_single, [Nb(2), 1]), 2);
                    tN_mean = mean(tN_single, 2);
                else
                    tY_weighted = zeros(Nb(2), 1);
                    tN_mean = 0;
                end

                for ne = 1:Nb(2)
                    Y_idx = (ne - 1) * numel(ss.VN) + global_fold_idx;
                    Y(Y_idx) = tY_weighted(ne);
                end
                N(global_fold_idx) = tN_mean;
            end
        end

        % Normalise by localizer weights, then average across cross-validation folds
        Y = reshape(Y ./ max(eps, repmat(N, [Nb(2), 1])), [size(N, 1), Nb(2)]);
        Y = ss.PV * Y;  % weighted average across folds per subject
        N = 1 ./ (ss.PV * (1 ./ max(eps, N)));  % harmonic mean of localizer weights
        sN = mean(N > 1e-4, 1);
    end

    if sN > 0
        % Estimate GLM: OLS or ReML depending on ss.estimation
        if strcmpi(ss.estimation, 'ols')
            iC = double(N > 1e-4);  % binary weights (handles missing data)
        else
            % ReML: estimate covariance from residuals
            n_weights = N;
            y = Y .* sqrt(n_weights(:, ones(1, Nb(2))));
            x = ss.X .* sqrt(n_weights(:, ones(1, Nb(1))));
            e = Y - ss.X * (pinv(x' * x) * (x' * y));
            [~, iC] = spm_ss_fls({e, n_weights});
        end
        % Whiten data and design matrix, then fit GLM
        y = Y .* iC(:, ones(1, Nb(2)));
        x = ss.X .* iC(:, ones(1, Nb(1)));
        [b, ee] = spm_ss_glm('estimate', x, y);
        Bplane(:, :, nroi) = b;
        Cplane(:, nroi) = iC;
        Eplane(:, :, nroi) = ee;
    end
    Pplane(nroi) = mean(N);
    Oplane(nroi) = sN;
    Zplane(:, :, nroi) = Y;
end
fprintf(1, '\n');

%% Write output volumes
% beta estimates: one volume per regressor x effect combination
nb=1;for nb1=1:Nb(1),for nb2=1:Nb(2),z=nan+zeros(size(frois{1}));for nroi=1:nrois,z(frois{1}==ss.VM_roiids(nroi))=Bplane(nb1,nb2,nroi);end; spm_write_vol(VB(nb),z);nb=nb+1;end;end
% overlap: proportion of subjects with data in each ROI
z=nan+zeros(size(frois{1}));for nroi=1:nrois,z(frois{1}==ss.VM_roiids(nroi))=Oplane(nroi);end; spm_write_vol(VO,z);

% save the ss structure with estimation results
ss.estimate=struct('BETA',VB,'OVERLAP',VO,'beta',Bplane,'rss',Eplane,'whitening',Cplane,'overlap',Oplane,'voxels',Nplane,'coverage',Pplane,'qa',Mplane,'y',Zplane);
save(fullfile(ss.swd,['SPM_ss',extname,'.mat']),'ss');
disp(['Analysis file saved: ',fullfile(ss.swd,['SPM_ss',extname,'.mat'])]);

%% Write CSV output files
% Build subject/fold index mapping for QA reporting
nidxs=zeros(ss.n,1);
nidxs2=zeros(2,size(ss.estimate.qa,1));
for nf=1:size(ss.estimate.qa,1),
    [nill,idxs]=max(ss.PV(:,nf));
    nidxs(idxs)=nidxs(idxs)+1;
    nidxs2(1,nf)=idxs;
    nidxs2(2,nf)=nidxs(idxs);
end

if 0, % note: this will be removed in the future; change to 1 to create old-format _data.csv file
    fname=['spm_ss',extname,'_data.csv'];
    fh=fopen(fullfile(ss.swd,fname),'wt');
    fprintf(fh,'Data\n');
    fprintf(fh,'ROI#,average ROI size,average localizer mask size,inter-subject overlap,');
    for ns=1:ss.n,for ne=1:Nb(2),fprintf(fh,'Subject#%d[%d]',ns,ne); if ne<Nb(2)||ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end
    for nroi=1:nrois,
        fprintf(fh,'%d,%d,%d,%f,',nroi,round(mean(ss.estimate.voxels(:,nroi))),round(mean(ss.estimate.voxels(:,nroi))*ss.estimate.coverage(nroi)),ss.estimate.overlap(nroi));
        for ns=1:ss.n, for ne=1:Nb(2),fprintf(fh,'%f',Zplane(ns,ne,nroi)); if ne<Nb(2)||ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end
    end
    fprintf(fh,'\nWeights\n');
    fprintf(fh,'ROI#,');
    for ns=1:ss.n,fprintf(fh,'Subject#%d',ns); if ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end;
    for nroi=1:nrois,
        fprintf(fh,'%d,',nroi);
        for ns=1:ss.n, fprintf(fh,'%f',Cplane(ns,nroi)); if ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end;
    end
    fprintf(fh,'\nquality control (localizer mask sizes)\n');
    fprintf(fh,'Subject#,Session/Partition#,filename'); for nroi=1:nrois,fprintf(fh,',ROI#%d',nroi);end;fprintf(fh,'\n');
    for nf=1:size(ss.estimate.qa,1),
        fprintf(fh,'%d,%d,%s',nidxs2(1,nf),nidxs2(2,nf),ss.PN{nf});
        for nroi=1:nrois,fprintf(fh,',%d',round(ss.estimate.qa(nf,nroi)));end
        fprintf(fh,'\n');
    end
    fclose(fh);
end

% New format CSV output files
if size(ss.EffectOfInterest_contrasts,2)==Nb(2), effect_names=ss.EffectOfInterest_contrasts;
else effect_names=arrayfun(@num2str,1:Nb(2),'uni',0);
end
fname=['spm_ss',extname,'_data.details.Weight.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Subject,Weight\n');
for nroi=1:nrois,for ns=1:ss.n, fprintf(fh,'%s,%s,%f\n',ss.VM_roinames{nroi},ss.PV_subjnames{ns},Cplane(ns,nroi)); end; end
fclose(fh);
fname=['spm_ss',extname,'_data.summaries.EffectSize.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Effect,MeanEffect,StdEffect,StderrEffect\n');
for nroi=1:nrois,for ne=1:Nb(2), fprintf(fh,'%s,%s,%f,%f,%f\n',ss.VM_roinames{nroi},effect_names{1,ne},mean(Zplane(:,ne,nroi)),std(Zplane(:,ne,nroi)),std(Zplane(:,ne,nroi))/sqrt(size(Zplane,1))); end; end
fclose(fh);
fname=['spm_ss',extname,'_data.details.ROIsize.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Subject,Session/Partition,LocalizerSize\n');
for nroi=1:nrois, for nf=1:size(ss.estimate.qa,1), fprintf(fh,'%s,%s,%d,%d\n',ss.VM_roinames{nroi},ss.PV_subjnames{nidxs2(1,nf)},nidxs2(2,nf),round(ss.estimate.qa(nf,nroi)));end; end
fclose(fh);
fname=['spm_ss',extname,'_data.summaries.ROIsize.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,ROISize,LocalizerSize,LocalizerIntersubjectOverlap\n');
for nroi=1:nrois, fprintf(fh,'%d,%d,%d,%f\n',ss.VM_roinames{nroi},round(mean(ss.estimate.voxels(:,nroi))),round(mean(ss.estimate.voxels(:,nroi))*ss.estimate.coverage(nroi)),ss.estimate.overlap(nroi)); end
fname=['spm_ss',extname,'_data.details.SourceFiles.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'Subject,Session/Partition,Localizer');
for ne=1:Nb(2), fprintf(fh,',%s',effect_names{1,ne}); end; fprintf(fh,'\n');
for nf=1:size(ss.estimate.qa,1),
    fprintf(fh,'%s,%d,%s',ss.PV_subjnames{nidxs2(1,nf)},nidxs2(2,nf),ss.PN{nf});
    for ne=1:size(ss.PY,2), fprintf(fh,',%s',ss.PY{nf,ne}); end; fprintf(fh,'\n');
end
fclose(fh);
fname=['spm_ss',extname,'_data.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Subject,Effect, LocalizerSize,EffectSize\n');
for nroi=1:nrois,for ns=1:ss.n, for ne=1:Nb(2), fprintf(fh,'%s,%s,%s, %d,%f\n',ss.VM_roinames{nroi},ss.PV_subjnames{ns},effect_names{1,ne}, round(ss.PV(ns,:)*ss.estimate.qa(:,nroi)),Zplane(ns,ne,nroi)); end; end; end
fclose(fh);

%% Generate QA plots (requires CONN toolbox)
try, spm_ss qacreate;
catch, fprintf('warning: unable to generate QA plots (possibly missing graphic display capabilities). Please use "spm_ss qa ''%s''" syntax to create these plots at a later time\n',ss.swd);
end

%% Evaluate defined contrasts
ss=spm_ss_contrast_ROI(ss);
cd(cwd);

end


function conn_savetextfile(tfilename,data,names,descrip)
% conn_savetextfile saves numeric data data to text file
% conn_savetextfile(tfilename,data [,names,descrip])
%

if nargin<4||isempty(descrip), descrip={}; end
if nargin<3||isempty(names), names={}; end
[nill,nill,tfileext]=fileparts(tfilename);
switch(tfileext)
    case '.mat'
        if ~isempty(names)&&~isempty(descrip), save(tfilename,'data','names','descrip');
        elseif ~isempty(names), save(tfilename,'data','names');
        else save(tfilename,'data');
        end
    otherwise,
        if strcmp(tfileext,'.txt'), names=regexprep(names,'\s','');
        else                        names=regexprep(names,'\,','');
        end
        fh=fopen(tfilename,'wt');
        for n1=1:numel(names),
            if isempty(names{n1}), names{n1}='-'; end
            fprintf(fh,'%s',names{n1});
            if n1<numel(names)&&strcmp(tfileext,'.csv'), fprintf(fh,',');
            elseif n1<numel(names), fprintf(fh,' ');
            else fprintf(fh,'\n');
            end
        end
        for n2=1:size(data,1),
            for n1=1:size(data,2),
                if iscell(data(n2,n1))&&ischar(data{n2,n1}), fprintf(fh,'%s',data{n2,n1});
                else fprintf(fh,'%s',mat2str(data(n2,n1)));
                end
                if n1<size(data,2)&&strcmp(tfileext,'.csv'), fprintf(fh,',');
                elseif n1<size(data,2), fprintf(fh,' ');
                else fprintf(fh,'\n');
                end
            end
        end
        fclose(fh);
end
end
