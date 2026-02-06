function flip_image_sagittal(mask_fname, output_fname)
% FLIP_MASK_SAGITTAL Flips a NIfTI mask along the sagittal plane
%
% Inputs:
%   mask_fname - Filename of input mask
%   output_fname - Filename for flipped output mask (optional)
%
% Example:
%   flip_mask_sagittal('/path/to/dir', 'mask.nii', 'flipped_mask.nii')

% Handle optional output filename
if nargin < 3
    [~, name, ext] = fileparts(mask_fname);
    output_fname = [name '_flipped' ext];
end

% Load and flip mask
mask_vol = spm_vol(mask_fname);
mask_data = spm_read_vols(mask_vol);
flipped_mask_data = flip(mask_data, 1);

% Create output volume and save
flipped_vol = mask_vol;
flipped_vol.fname = output_fname;
spm_write_vol(flipped_vol, flipped_mask_data);

end