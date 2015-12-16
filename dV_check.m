%close all
clear all

addpath C:\Users\erikakl\Documents\MATLAB\spm12\
close all
warning off all
spm('Defaults','fmri')
global defaults
defaults.analyze.flip = 1;

folder_subjects = strcat(pwd,'/Subjects_nifti');
subject_folders = dir(folder_subjects);
subject_folders = subject_folders(3:end);

velocity_name = 'ccj.nii';
%mask_name = '51_mask.nii';
%mask_name = '51_mask_no1.nii';
mask_name = 'ccj_mask.nii';
dicom_name = 'hdr.dcm';

total_delta_P = [];
teller = 1;

%a = load('sub18_meanvel_mm.txt');



Q_all = zeros(9,3);
for p = 22
    
    teller = teller+1;
    disp(teller)
    current_folder = strcat(folder_subjects,'/',subject_folders(p).name);
    file_status_ccj = dir(fullfile(current_folder,velocity_name));                           
    
    %if  isempty(file_status_ccj)==0                     
     %   k
        fullfile(current_folder,velocity_name);    
        velocity_volume = spm_vol(fullfile(current_folder,velocity_name));
        mask_volume = spm_vol(fullfile(current_folder,mask_name));

%velocity_image_file = fullfile('./Subject_016/ccj.nii');
%mask_image_file = fullfile('./Subject_016/ccj_mask.nii');
%current_folder = './Subject_016';
%dicom_name = '2.dcm';

density = 1.0007; % density for CSF is 1.0007 g/(cm^3), (this value comes from Alperin)
visc = 1.1*0.01;%mu for CSF is about 1.1 cP ("centipoise"), where 1 cP = 0.01 g/(cm*s) (this value comes from Alperin)


%velocity_volume = spm_vol(velocity_image_file);
%mask_volume = spm_vol(mask_image_file);


velocity_image = zeros(velocity_volume.dim(1),velocity_volume.dim(2),velocity_volume.dim(3)); 
mask_image = flipdim(spm_slice_vol(mask_volume,spm_matrix([0 0 1]),mask_volume.dim(1:2),1),2);
mask_image = im2bw(mask_image,1e-6);
%NB!!!!!!!!!!!!!!!!!!! Kun for 47 og 51!!!!!
mask_image = fliplr(mask_image);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = dicominfo(fullfile(current_folder,dicom_name));
final_phases = velocity_volume.dim(3);%time
dx = info.PixelSpacing(1)*0.1;%mm to cm
dy = info.PixelSpacing(2)*0.1;%mm to cm
dt = 60/info.HeartRate/final_phases;  %Heart rate of object, final phases: t 
V_enc = info.Private_2001_101a(3);%Max value possible
%disp(V_enc) 
alf = 1.1; %why? Used in filter
thresh1 = 2.5; %why? Used in filter

idx = find(mask_image);
[xx,yy] = meshgrid(1:info.Width);

figure;
plot(xx(idx),yy(idx),'c.');

%spm_slice_vol returns a section through a memory mapped image volume on
%disk. This section is the transverse slize at z=0 after linear
%transformation according to matrix A
for n = 1:velocity_volume.dim(3)
    velocity_image(:,:,n) = spm_slice_vol(velocity_volume,spm_matrix([0 0 n]),velocity_volume.dim(1:2),1);
end 

numpix = length(find(mask_image>0)); 

%Size of ROI: numpix*dx*dy
disp(numpix*dx*dy);

svv = size(velocity_image);
velocity_integrals = zeros(1,svv(3));
velocity_integrals_pos = zeros(1,svv(3));
velocity_integrals_neg = zeros(1,svv(3));
Q = zeros(1,svv(3));
time = 0:dt:(svv(3)-1)*dt;

vk = zeros(svv(3),1);
intk = zeros(svv(3),1);
maxk = zeros(svv(3),1);
mink = zeros(svv(3),1);
antzero = zeros(svv(3),1);
peakvel = zeros(svv(3),2);


for k=1:svv(3)
    tmp = double(velocity_image(:,:,k));
    intk(k) = sum(tmp(idx))*dx*dy;
    maxk(k) = max(tmp(idx));
    mink(k) = min(tmp(idx));
    vk(k) = sum(tmp(idx));
    sc = size(find(tmp(idx)==0));
    antzero(k) = sc(1,1)/numpix;
    peakvel(k,1) = max(max(tmp(idx)));
    peakvel(k,2) = min(min(tmp(idx)));
    
    
    
end
% 
ci = 23;
figure;imagesc(velocity_image(:,:,ci))
colormap(gray)
colorbar
hold on
plot(xx(idx),yy(idx),'r.');
title(subject_folders(p).name)
% 

figure;imagesc(velocity_image(:,:,ci))
colormap(gray)
colorbar
title(subject_folders(p).name)

%figure;imagesc(velocity_image(:,:,ci))
%colormap(gray)
%colorbar
%for l = 10
%figure;
%hold on
figure;imagesc(velocity_image(:,:,ci).*double(mask_image));
%axis([100 150 100 180])
colorbar
title(subject_folders(p).name)
%colorscale(gray)
% %%%end

% figure
% plot(tmp(idx),'.')

% vk1 = vk;
vk = vk/(length(idx));

vkcm = zeros(svv(3),1);
maxkcm = vkcm;
minkcm = vkcm;
peakvelcm = zeros(svv(3),2);
%convert pix-values to cm/s
%4.094=10cm/s, 2.047=0cm/s, 0=-10cm/s
for i = 1:svv(3)
    vkcm(i) = 20*(vk(i)-2.047)/4.094;
    maxkcm(i) = 20*(maxk(i)-2.047)/4.094;
    minkcm(i) = 20*(mink(i)-2.047)/4.094;
    peakvelcm(i,1) = 20*(peakvel(i,1)-2.047)/4.094;
    peakvelcm(i,2) = 20*(peakvel(i,2)-2.047)/4.094;
end


%vk2 = vk/(length(idx)) + 0.4*antzero; %0.5*antzero*(vk/numpix);%vk/(numpix-150) + .2;%((vk/length(idx))+ 0.4*(vk./numpix)) + .2;%;


figure
plot(time,vkcm,'bo-')
title(subject_folders(p).name)
axis([0 1 -5 5])
%vk2 = vk;

% figure
% plot(vk)
% 
% x = 1:length(a);
% x2 = 2:svv(3)+1;
% % 

% figure
% plot(x,a,'o-')
% hold on
% plot(x2-1,(vk*10),'ro-')  
% title([subject_folders(p).name, 'mean velocity'])
% ylabel('mean vel. [cm/s]')
% legend('Master','calc.')
% 
% figure
% plot(x,a,'o-')
% hold on
% plot(x2-1,(vk2*10)+0,'ro-',x2-1,vk*10,'mo-')
% title([subject_folders(p).name,'mean velocity, adjusted'])
% ylabel('mean vel. [cm/s]')
% legend('Master','calc. adj','calc')
% 
% abspeak = zeros(1,svv(3));
% Tmp = 0;
% 
% for i = 1:length(abspeak)
%     Tmp = abs(peakvel(i,1));
%     if abs(peakvel(i,2)) > Tmp
%         Tmp = abs(peakvel(i,2));
%     end
%     abspeak(i) = Tmp;
% end

% figure
% plot(x2-1,peakvel(:,1),'bo-',x2-1,peakvel(:,2),'ro-',x2-1,abspeak,'ko-')
% legend('+','-','abs')    
% 
% pv = load('sub18_peakvel_mm.txt');
% 
% figure
% plot(x,pv,'bo-',x2-1,abspeak*10,'ro-')
% % legend('Geir','Erika')
% ylabel('mm/s')
% 
% % 
% A = load('sub18_minvel_mm.txt');
% B = load('sub18_maxvel_mm.txt');
% 
% figure
% plot(x,A,'bo-',x2,mink*10,'ro-')
% legend('Geir','Erika');
% title('Min-vel')
% ylabel('mm/s')
% 
% figure
% plot(x,B,'bo-',x2,maxk*10,'ro-')
% legend('Geir','Erika');
% title('Max-vel')
% ylabel('mm/s')
% 
% figure
% plot(x2,mink*10,'c',x2,maxk*10,'m')
% % figure
% % plot(a,'co-')
% % hold on
% %plot(mink*10,'mo-')
% 
% a = load('flux_sub18_mm.txt');
% figure
% plot(x,a,'b.-',x2,(intk*1000)+00,'r.-')
% legend('Geir','Erika')
% title('flux')
% ylabel('mm^3/s')

%figure
%plot(x2, intk*1000)





end