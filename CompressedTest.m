close all

folder_subjects = strcat(pwd,'/Subjects_nifti');
subject_folders = dir(folder_subjects);
subject_folders = subject_folders(3:end);

velocity_name = 'ccj.nii';
mask_name = 'ccj_mask.nii';
dicom_name = 'hdr.dcm';

%Loop for all subjects.
for p = 1
    
current_folder = strcat(folder_subjects,'/',subject_folders(p).name);
file_status_ccj = dir(fullfile(current_folder,velocity_name));                           
velocity_volume = spm_vol(fullfile(current_folder,velocity_name));
mask_volume = spm_vol(fullfile(current_folder,mask_name));

density = 1.0007; % density for CSF is 1.0007 g/(cm^3), (this value comes from Alperin)
visc = 1.1*0.01;%mu for CSF is about 1.1 cP ("centipoise"), where 1 cP = 0.01 g/(cm*s) (this value comes from Alperin)

velocity_image = zeros(velocity_volume.dim(1),velocity_volume.dim(2),velocity_volume.dim(3)); 
mask_image = flipdim(spm_slice_vol(mask_volume,spm_matrix([0 0 1]),mask_volume.dim(1:2),1),2);
%Convert ROI to black and white.
mask_image = im2bw(mask_image,1e-6);


info = dicominfo(fullfile(current_folder,dicom_name));
final_phases = velocity_volume.dim(3);%time
dx = info.PixelSpacing(1)*0.1;%mm to cm
dy = info.PixelSpacing(2)*0.1;%mm to cm
dt = 60/info.HeartRate/final_phases;  %Heart rate of object, final phases: t 
V_enc = info.Private_2001_101a(3);%Max value possible, if exceeded: aliasing

%idx corresponds to the ROI-position. Where pixels != 0 in mask_image
%xx and yy are used for visaulizing ROI.
idx = find(mask_image);
[xx,yy] = meshgrid(1:info.Width);

%spm_slice_vol returns a section through a memory mapped image volume on
%disk. This section is the transverse slize at z=0 after linear
%transformation according to matrix A
for n = 1:velocity_volume.dim(3)
    velocity_image(:,:,n) = spm_slice_vol(velocity_volume,spm_matrix([0 0 n]),velocity_volume.dim(1:2),1);
end 


end