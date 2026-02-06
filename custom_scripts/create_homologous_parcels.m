function create_homologous_parcels(input_fname, output_fname, smooth_kernel, overlap_threshold)
% CREATE_HOMOLOGOUS_PARCELS Creates bilateral, homologous ROIs based on an overlap map
%
% Inputs:
%   input_fname - Filename of input overlap map (.nii)
%   output_fname - Filename for output ROIs (.nii) (optional)
%   smooth_kernel - FWHM in mm for smoothing [x y z] (optional, default [6 6 6])
%   overlap_threshold - Threshold for including voxels (optional, default 0.2)
%
% Example:
%   create_homologous_parcels('Overlap.nii', 'fROIs_homologous.nii', [8 8 8], 0.1)

% Handle optional arguments
if nargin < 2
    [path, name, ~] = fileparts(input_fname);
    output_fname = fullfile(path, [name '_homologous_rois.nii']);
end

if nargin < 3
    smooth_kernel = [6 6 6];
end

if nargin < 4
    overlap_threshold = 0.2;
end

% 1. Read in the original overlap map
fprintf('Loading original overlap map: %s\n', input_fname);
overlap_vol = spm_vol(input_fname);
overlap_data = spm_read_vols(overlap_vol);

% 2. Create a temporary file for the flipped map
[path, name, ext] = fileparts(input_fname);
flipped_fname = fullfile(path, [name '_flipped' ext]);

% 3. Flip the overlap map across the sagittal plane
fprintf('Flipping overlap map across sagittal plane...\n');
flip_image_sagittal(input_fname, flipped_fname);

% 4. Load the flipped map
flipped_vol = spm_vol(flipped_fname);
flipped_data = spm_read_vols(flipped_vol);

% 5. Create homologous map by averaging original and flipped maps
fprintf('Creating homologous map by averaging original and flipped maps...\n');
homologous_data = (overlap_data + flipped_data) / 2;

% 6. Save the homologous map to a temporary file
homologous_fname = fullfile(path, [name '_homologous' ext]);
homologous_vol = overlap_vol;
homologous_vol.fname = homologous_fname;
spm_write_vol(homologous_vol, homologous_data);

% 7. Smooth the homologous map
fprintf('Smoothing homologous map with FWHM = [%d %d %d]...\n', ...
    smooth_kernel(1), smooth_kernel(2), smooth_kernel(3));
smooth_fname = fullfile(path, ['s' name '_homologous' ext]);
spm_smooth(homologous_fname, smooth_fname, smooth_kernel);

% 8. Read in the smoothed map
smooth_vol = spm_vol(smooth_fname);
smooth_data = spm_read_vols(smooth_vol);

% 9. Find voxels that meet overlap threshold
fprintf('Finding voxels above threshold (%.2f)...\n', overlap_threshold);
voxels_above_threshold = find(smooth_data >= overlap_threshold);

% 10. Run watershed on negative of smoothed map to create ROIs
fprintf('Running watershed segmentation...\n');
segmented_regions = spm_ss_watershed(-smooth_data, voxels_above_threshold);

% 11. Save segmented regions
fprintf('Saving segmented ROIs to: %s\n', output_fname);
output_vol = smooth_vol; % Clone header info
output_vol.fname = output_fname;
output_vol.dt = [spm_type('int16') spm_platform('bigend')]; % Change datatype to int16 for ROI labels
spm_write_vol(output_vol, segmented_regions);

% 12. Clean up temporary files
fprintf('Cleaning up temporary files...\n');
delete(flipped_fname);
delete(homologous_fname);
delete(smooth_fname);

fprintf('Homologous ROIs created successfully.\n');
fprintf('Number of ROIs created: %d\n', max(segmented_regions(:)));
end
