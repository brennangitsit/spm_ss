function reindexed_nii = reindex_homologous_parcels(input_nii, parcel_data, varargin)
% REINDEX_HOMOLOGOUS_PARCELS - Reindex parcels with homologous regions
%
% This function zeros out the sagittal midline of a NIfTI image and reindexes
% parcels according to a data table with hemisphere information, preserving the
% exact order of parcels as provided in the input table.
%
% Inputs:
%   input_nii   - Path to input NIfTI file or loaded NIfTI structure
%   parcel_data - Table or structure with fields:
%                 'parcel' - original parcel numbers
%                 'hemi'   - hemisphere ('LH', 'RH', or 'BI')
%                 'label'  - parcel labels
%   midline     - (Optional) Custom midline specification:
%                 - Single number: Index of the midline slice to remove
%                 - Two numbers [a b]: Use indices a and b as the boundaries
%                   for left and right hemispheres (for even dimensions)
%                 - If not provided, automatically uses round(dimensions(1)/2)
%
% Output:
%   reindexed_nii - Path to the saved reindexed NIfTI file
%
% Examples:
%   reindexed_nii = reindex_homologous_parcels('parcels.nii', parcel_table);
%   reindexed_nii = reindex_homologous_parcels('parcels.nii', parcel_table, 27);
%   reindexed_nii = reindex_homologous_parcels('parcels.nii', parcel_table, [26 27]);

% Initialize MarsBaR for NIfTI operations (if it's not already initialized)
if ~exist('marsbar', 'file') || ~isappdata(0, 'mars_vars')
    try
        marsbar('on');
    catch
        warning('MarsBaR not found. Some functions might not work properly.');
    end
end

% Load the input NIfTI image
if ischar(input_nii) || isstring(input_nii)
    vol = spm_vol(input_nii);
    fprintf('Loaded NIfTI file: %s\n', input_nii);
else
    vol = input_nii;
    fprintf('Using provided NIfTI structure\n');
end
img = spm_read_vols(vol);

% Store the original dimensions
orig_dims = size(img);
fprintf('Image dimensions: %dx%dx%d\n', orig_dims);

% Convert parcel_data to table if it's not already
if ~istable(parcel_data)
    parcel_data = struct2table(parcel_data);
end

% Make sure required fields exist
required_fields = {'parcel', 'hemi', 'label'};
for i = 1:length(required_fields)
    if ~ismember(required_fields{i}, parcel_data.Properties.VariableNames)
        error('Parcel data must contain field: %s', required_fields{i});
    end
end

% Display the input table for verification
fprintf('\nInput parcel table (%d parcels):\n', height(parcel_data));
disp(parcel_data);

% Check what values actually exist in the image
unique_values = unique(img(:));
unique_values = unique_values(unique_values > 0);  % Exclude zero
fprintf('Values found in the image: %s\n', mat2str(unique_values));

% Check that all parcels exist in the image
fprintf('\nValidating parcels...\n');
unique_img_parcels = unique(img(:));
unique_img_parcels = unique_img_parcels(unique_img_parcels > 0); % Exclude background (0)
unique_table_parcels = unique(parcel_data.parcel);

% Check if all parcels in the table exist in the image
missing_parcels = setdiff(unique_table_parcels, unique_img_parcels);
if ~isempty(missing_parcels)
    warning('The following parcels in the table are not present in the image: %s', ...
        mat2str(missing_parcels));
end

% Check if there are parcels in the image not listed in the table
unlisted_parcels = setdiff(unique_img_parcels, unique_table_parcels);
if ~isempty(unlisted_parcels)
    warning('The following parcels are in the image but not listed in the table: %s', ...
        mat2str(unlisted_parcels));
end

% Handle the midline parameter
if nargin >= 3 && ~isempty(varargin{1})
    midline_param = varargin{1};
    
    % Check if midline parameter is valid
    if numel(midline_param) > 2
        error('The midline parameter must be either a single value or a vector of two values.');
    end
    
    % Check if provided indices are within image dimensions
    if any(midline_param < 1) || any(midline_param > orig_dims(1))
        error('Midline indices must be within the image dimensions (1 to %d).', orig_dims(1));
    end
    
    if numel(midline_param) == 1
        % Single midline value - remove this slice
        midline_idx = midline_param;
        lh_boundary = midline_idx - 1;  % Last index of LH
        rh_boundary = midline_idx + 1;  % First index of RH
        has_midline_slice = true;
        fprintf('Using custom midline at x = %d\n', midline_idx);
    else
        % Two values - boundaries between hemispheres
        lh_boundary = midline_param(1);
        rh_boundary = midline_param(2);
        has_midline_slice = false;
        fprintf('Using custom hemisphere boundaries: LH ends at x = %d, RH starts at x = %d\n', ...
            lh_boundary, rh_boundary);
    end
else
    % Auto determine midline based on dimensions
    midline_idx = round(orig_dims(1)/2);
    lh_boundary = midline_idx - 1;
    rh_boundary = midline_idx + 1;
    has_midline_slice = true;
    fprintf('Using auto-detected midline at x = %d\n', midline_idx);
end

% Create left and right hemisphere masks
lh_mask = false(orig_dims);
rh_mask = false(orig_dims);
midline_mask = false(orig_dims);

lh_mask(1:lh_boundary,:,:) = true;      % LH is at lower X indices
rh_mask(rh_boundary:end,:,:) = true;    % RH is at higher X indices

if has_midline_slice
    midline_mask(midline_idx,:,:) = true;  % The midline slice
end

% Save masks for visualization if needed
mask_vol = vol;
mask_vol.fname = 'lh_mask.nii';
spm_write_vol(mask_vol, single(lh_mask));
mask_vol.fname = 'rh_mask.nii';
spm_write_vol(mask_vol, single(rh_mask));
if has_midline_slice
    mask_vol.fname = 'midline_mask.nii';
    spm_write_vol(mask_vol, single(midline_mask));
end

% Create a new empty volume for the reindexed parcels
img_reindexed = zeros(orig_dims);

% Initialize a table for the new mapping
new_mapping = table();
new_idx_counter = 1;

fprintf('\nProcessing parcels in the order provided...\n');
% Process the parcels in the exact order they appear in the input table
for i = 1:height(parcel_data)
    orig_idx = parcel_data.parcel(i);
    label = parcel_data.label{i};
    hemi = parcel_data.hemi{i};
    parcel_mask = (img == orig_idx);
    
    % Hemisphere-specific masks
    switch hemi
        case 'BI'
            hemispheres = {'LH', lh_mask; 'RH', rh_mask};
        otherwise
            hemispheres = {hemi, strcmp(hemi, 'LH') * lh_mask + strcmp(hemi, 'RH') * rh_mask};
    end

    for h = 1:size(hemispheres,1)
        hemi_tag = hemispheres{h,1};
        mask = hemispheres{h,2};
        hemi_masked = parcel_mask & mask;
        
        if any(hemi_masked(:))
            row_idx = height(new_mapping) + 1;
            new_mapping.parcel(row_idx) = new_idx_counter;
            new_mapping.hemi_label{row_idx} = [hemi_tag '_' label];
            new_mapping.hemi{row_idx} = hemi_tag;
            new_mapping.label{row_idx} = label;
            new_mapping.original_parcel(row_idx) = orig_idx;
            img_reindexed(hemi_masked) = new_idx_counter;
            new_idx_counter = new_idx_counter + 1;
            fprintf('Parcel %d (%s, %s): assigned index %d\n', orig_idx, hemi_tag, label, new_idx_counter - 1);
        else
            fprintf('Parcel %d (%s, %s): no voxels found\n', orig_idx, hemi_tag, label);
        end
    end
end

% Save the reindexed image
[filepath, name, ext] = fileparts(vol.fname);
output_fname = fullfile(filepath, [name '_reindexed' ext]);
vol_out = vol;
vol_out.fname = output_fname;
spm_write_vol(vol_out, img_reindexed);

% Display the new mapping
fprintf('\nNew parcel mapping (%d parcels):\n', height(new_mapping));
disp(new_mapping);

% Save the new mapping to a CSV file
mapping_file = fullfile(filepath, [name '_new_mapping.csv']);
writetable(new_mapping, mapping_file);

% Return the output file path
reindexed_nii = output_fname;

fprintf('\nReindexed image saved as: %s\n', output_fname);
fprintf('New parcel mapping saved as: %s\n', mapping_file);
end
