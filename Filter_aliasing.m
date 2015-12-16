close all

addpath C:\Users\erikakl\Documents\MATLAB\spm12\


folder_subjects = strcat(pwd,'/Subjects_nifti');
subject_folders = dir(folder_subjects);
subject_folders = subject_folders(3:end);

velocity_name = 'ccj.nii';
mask_name = 'ccj_mask.nii';
dicom_name = 'hdr.dcm';

Q_all = zeros(9,3);
for p = 13
    
current_folder = strcat(folder_subjects,'/',subject_folders(p).name);
file_status_ccj = dir(fullfile(current_folder,velocity_name));                           
velocity_volume = spm_vol(fullfile(current_folder,velocity_name));
mask_volume = spm_vol(fullfile(current_folder,mask_name));

density = 1.0007; % density for CSF is 1.0007 g/(cm^3), (this value comes from Alperin)
visc = 1.1*0.01;%mu for CSF is about 1.1 cP ("centipoise"), where 1 cP = 0.01 g/(cm*s) (this value comes from Alperin)

velocity_image = zeros(velocity_volume.dim(1),velocity_volume.dim(2),velocity_volume.dim(3)); 
mask_image = flipdim(spm_slice_vol(mask_volume,spm_matrix([0 0 1]),mask_volume.dim(1:2),1),2);
mask_image = im2bw(mask_image,1e-6);
%mask_image = fliplr(mask_image);

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
    tmp = ((2*V_enc)/4.094)*tmp-V_enc;
    balance(k,1) = size(find(tmp(idx)>0),1);
    balance(k,2) = size(find(tmp(idx)<0),1);
    balance(k,3) = ((balance(k,1)+balance(k,2))-numpix);
    balance(k,4) = numpix;
    %figure; imagesc(tmp*mask_image); colorbar
end
points = 1:2;%1:svv(3);

%aviobj = VideoWriter('mymovie.avi'); 
xylim = [145 185 120 165];


%[120 160 100 160];%2, 4
%[130 165 115 160];%9
%[105 140 95 155];%10
%[130 160 115 160];%11
%[115 145 105 150];%12
%[100 150 105 165]%14
%[110 150 105 155];%16
%[125 160 90 140]%17
%[110 150 110 160]%18
%[120 165 95 150]%19
%[135 175 105 155]%20
%[145 185 120 165]%22
%[120 155 110 155];%23
%[150 185 100 155];%25
%[110 145 100 150];%28
%[130 170 100 155]%29
%[155 195 105 155]%33
%[110 150 105 160 ]%41
%[130 180 120 180]%46
%[140 180 115 165];%47


%aviobj1 = avifile('CSF_Subject22.avi','fps',2); 

std_vk = zeros(svv(3),1);
TMP = zeros(svv(1),svv(2),svv(3));
for k=points;%:svv(3)
    tmp = double(velocity_image(:,:,k));
    tmp = ((2*V_enc)/4.094)*tmp-V_enc;
    TMP(:,:,k) = (tmp.*double(mask_image));
%     h = imagesc(TMP(:,:,k).*double(mask_image));colorbar;colormap('default');%colormap(gray)
%     title([num2str(time(k),'%3d'),' s   Subject ',subject_folders(p).name])
%     axis(xylim);caxis([-1 1]*6)
%     set(gca, 'CLim',[-1 1]*6);
%     set(h,'EraseMode','xor');
%     frame = getframe(gca);
%     aviobj1 = addframe(aviobj1,frame);

    figure;imagesc(TMP(:,:,k).*double(mask_image));colorbar;colormap('default');%colormap(gray)
    title([num2str(k),' frame, Raw  Subject ',subject_folders(p).name])
    axis(xylim);caxis([-1 1]*10)
%     
end
%aviobj1 = close(aviobj1);

th_min = -2;%-3.5;
th_max = 2;%3.5;
TMP2 = TMP;

%aviobj = avifile('CSF_Subject22_Filtered.avi','fps',2); 

for k = points;
    for i = 1:svv(2)
        for j = 1:svv(2)
            if abs(peakvel(k,1))> 5.5 && abs(TMP(i,j,k)) ~=0 || abs(peakvel(k,2))> 5.5 && abs(TMP(i,j,k)) ~=0 
                
                if vk(k)>0 && TMP(i,j,k)<th_min
                    TMP2(i,j,k) = TMP2(i,j,k) + 2*V_enc; 
                end
                if vk(k)<0 && TMP(i,j,k)>th_max
                    TMP2(i,j,k) = TMP2(i,j,k) - 2*V_enc; 
                end
            end
            
        end
    end
    figure;imagesc(TMP2(:,:,k).*double(mask_image));colorbar;colormap('default');%colormap(gray)
    title([num2str(k),' frame, filtered   Subject ',subject_folders(p).name])
    axis(xylim);caxis([-1 1]*10)
    
%     h = imagesc(TMP2(:,:,k).*double(mask_image));colorbar;colormap('default');%colormap(gray)
%     title([num2str(k),' s   Subject ',subject_folders(p).name])
%     axis(xylim);caxis([-1 1]*6)
%     set(gca, 'CLim',[-1 1]*6);
%     set(h,'EraseMode','xor');
%     frame = getframe(gca);
%     aviobj = addframe(aviobj,frame);
end

%aviobj = close(aviobj);


% figure
% imagesc(TMP2(:,:,k));colorbar
% axis(xylim);caxis([-1 1]*6)
% title('Filter Test')
% pos = balance(:,1)./balance(:,4);
% figure
% plot(time(points),pos(points),'mp-',time,balance(:,2)./balance(:,4),'kp-')
% 
% 
% 

% % 
figure
plot(1:length(vk),vk,'ko-',1:length(vk),peakvel(:,1),'r*-',1:length(vk),peakvel(:,2),'b*-')
legend('Meanvel','MaxPos','MaxNeg')
axis([0 40 -6 6])
% 
% figure
% plot(1:length(vk),vk,'ko-',1:length(vk),peakvel(:,1),'r*-',1:length(vk),peakvel(:,2),'b*-',1:length(vk),peakvel(:,1)-std_vk,'rp-',1:length(vk),peakvel(:,2)+std_vk,'bp-')
% legend('Meanvel','MaxPos','MaxNeg','Max-std','Min-std')
% axis([0 40 -6 6])


% figure
% plot(1:length(vk),vk,'ko-',1:length(vk),vk+ vk*2.5,'r*-',1:length(vk),vk-vk*2.5,'b*-')
% legend('Meanvel','+50%','-50%')
% axis([0 40 -8 8])
% hold on
% plot(time,b(1:end-1)/10,'r*-');
% title(subject_folders(p).name)
% legend('Calc','Master')
% 
% figure
% plot(time,intk,'ro-');
% hold on
% plot(time,a(1:end-1)/1000,'r*-')%,time,b(1:end-1)/10,'r')
% title('Flux');
% 



% figure
% plot(time,vk,'bo-',time,intk,'ro-')
% title(subject_folders(p).name)
% axis([0 1 -7 5])
% legend('MeanVel','Flux')
end
