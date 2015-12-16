close all

addpath C:\Users\erikakl\Documents\MATLAB\spm12\


folder_subjects = strcat(pwd,'/Subjects_nifti');
subject_folders = dir(folder_subjects);
subject_folders = subject_folders(3:end);

velocity_name = 'ccj.nii';
mask_name = 'ccj_mask.nii';
dicom_name = 'hdr.dcm';

a = load('flux_sub18_mm.txt');
b = load('sub18_meanvel_mm.txt');

Q_all = zeros(9,3);
for p = 10
    
current_folder = strcat(folder_subjects,'/',subject_folders(p).name);
file_status_ccj = dir(fullfile(current_folder,velocity_name));                           
velocity_volume = spm_vol(fullfile(current_folder,velocity_name));
mask_volume = spm_vol(fullfile(current_folder,mask_name));

density = 1.0007; % density for CSF is 1.0007 g/(cm^3), (this value comes from Alperin)
visc = 1.1*0.01;%mu for CSF is about 1.1 cP ("centipoise"), where 1 cP = 0.01 g/(cm*s) (this value comes from Alperin)

velocity_image = zeros(velocity_volume.dim(1),velocity_volume.dim(2),velocity_volume.dim(3)); 
mask_image = flipdim(spm_slice_vol(mask_volume,spm_matrix([0 0 1]),mask_volume.dim(1:2),1),2);

info = dicominfo(fullfile(current_folder,dicom_name));
final_phases = velocity_volume.dim(3);%time
dx = info.PixelSpacing(1)*0.1;%mm to cm
dy = info.PixelSpacing(2)*0.1;%mm to cm
dt = 60/info.HeartRate/final_phases;  %Heart rate of object, final phases: t 
V_enc = info.Private_2001_101a(3);%Max value possible

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

numpix = length(find(mask_image>0)); 
%Size of ROI cm^2: numpix*dx*dy
disp(numpix*dx*dy);

%Store mean velocity (vk), volumetric flow rate (intk), number of positive and negative
%velocities within each time frame, min and max velcity within each time
%frame (peakvel)

svv = size(velocity_image);
time = 0:dt:(svv(3)-1)*dt;
velocity_integrals = zeros(1,svv(3));
Q = zeros(1,svv(3));
vk = zeros(svv(3),1);
intk = zeros(svv(3),1);
maxk = zeros(svv(3),1);
mink = zeros(svv(3),1);
peakvel = zeros(svv(3),2);
pos = zeros(svv(3),1);

%convert pix-values to cm/s
%4.094=6cm/s, 2.047=0cm/s, 0=-6cm/s
for k=1:svv(3)
    tmp = double(velocity_image(:,:,k));
    tmp = ((2*V_enc)/4.094)*tmp-V_enc;
    intk(k) = sum(tmp(idx))*dx*dy;
    vk(k) = sum(tmp(idx));
    peakvel(k,1) = max(max(tmp(idx)));
    peakvel(k,2) = min(min(tmp(idx)));
end
%Compute mean velocity, average over number of pixels
vk = vk/(length(idx));


balance = zeros(svv(3),4);
%Check percentange of + and - vel within each time frame. For fun, check
%flux also.
for k=1:svv(3)
    tmp = double(velocity_image(:,:,k));
    tmp = 20*(tmp-2.047)/4.094;
    balance(k,1) = size(find(tmp(idx)>0),1);
    balance(k,2) = size(find(tmp(idx)<0),1);
    balance(k,3) = ((balance(k,1)+balance(k,2))-numpix);
    balance(k,4) = numpix;
end



figure
plot(time,vk,'bo-',time,intk,'ro-')
title(subject_folders(p).name)
axis([0 1 -7 5])
legend('MeanVel','Flux')

figure
plot(time,vk,'bo-')
hold on
plot(time,b(1:end-1)/10,'r*-');
title(subject_folders(p).name)
legend('Calc','Master')

figure
plot(time,intk,'ro-');
hold on
plot(time,a(1:end-1)/1000,'r*-')%,time,b(1:end-1)/10,'r')
title('Flux');






end