close all
clear all

addpath C:\Users\erikakl\Documents\MATLAB\spm12\

folder_subjects = strcat(pwd,'/Subjects_nifti');
subject_folders = dir(folder_subjects);
subject_folders = subject_folders(3:end);

velocity_name = 'ccj.nii';
mask_name = 'ccj_mask.nii';
dicom_name = 'hdr.dcm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some of the subjects need special care:
% Remove noise pixels for subject 9, 10
% Mask Fliplr subject 47, 51
%Subject 4 (p=2) (143,117)
%Subject 9 (p=3), remove pix (138,135)
%Subject 10 (p=4), remove pix (106,125),(107,126)
%Subject 10 (p=4), remove pix (116,136)
%Loop for all subjects.
for p = 7
    
current_folder = strcat(folder_subjects,'/',subject_folders(p).name);
file_status_ccj = dir(fullfile(current_folder,velocity_name));                           
velocity_volume = spm_vol(fullfile(current_folder,velocity_name));
mask_volume = spm_vol(fullfile(current_folder,mask_name));

density = 1.0007; % density for CSF is 1.0007 g/(cm^3), (this value comes from Alperin)
visc = 1.1*0.01;%mu for CSF is about 1.1 cP ("centipoise"), where 1 cP = 0.01 g/(cm*s) (this value comes from Alperin)

N = 1;
velocity_image = zeros(velocity_volume.dim(1),velocity_volume.dim(2),velocity_volume.dim(3)); 
mask_image = flipdim(spm_slice_vol(mask_volume,spm_matrix([0 0 N]),mask_volume.dim(1:2),1),2);
%NB! Subj9 (p=3), (135,138)=0: p=4 (125,106),(126,107)=0
%mask_image(135,138) = 0;
%mask_image(106,125) =0; mask_image(107,126)=0;
%mask_Image(116,136) = 0; %p=1
%mask_image(122:129,122) = 0; %p=8
%Convert ROI to black and white.
mask_image = im2bw(mask_image,1e-6);
figure;imagesc(mask_image)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NB!!!!! Fliplr for 47 and 51
% Also for healthy subjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mask_image = fliplr(mask_image);

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

numpix = length(find(mask_image>0)); 
%Size of ROI cm^2: numpix*dx*dy
disp(numpix*dx*dy);

%Declare space for flow rate (Q), mean velocity (vk) and peak velocity
%(peakvel). 

svv = size(velocity_image);
time = 0:dt:(svv(3)-1)*dt;
Q = zeros(1,svv(3));
vk = zeros(svv(3),1);
peakvel = zeros(svv(3),2);

points =1:svv(3);

%Axis for figures.
xylim = [100 185 90 180];

TMP = zeros(svv(1),svv(2),svv(3));
for k=1:svv(3)
    tmp = double(velocity_image(:,:,k));
    tmp = ((2*V_enc)/4.094)*tmp-V_enc;
    TMP(:,:,k) = (tmp.*double(mask_image));     
    vk(k) = sum(tmp(idx));
    peakvel(k,1) = max(max(tmp(idx)));
    peakvel(k,2) = min(min(tmp(idx)));
end
vk = vk/(length(idx));

%Remove edge-pixels that contain noise.
if p == 2
    TMP(117,143,:) = 0;
end
if p ==3
    TMP(138,135,:) = 0;
end
if p==4
    TMP(106,125,:) = 0;
    TMP(107,126,:) = 0;
    TMP(116,136,:) = 0;
end

%Plot for safety check.
for k = 1%:svv(3)
 figure;imagesc(TMP(:,:,k).*double(mask_image));colorbar;colormap('jet');%colormap(gray)
    title([num2str(k),' frame, Raw  Subject ',subject_folders(p).name])
    axis(xylim);caxis([-1 1]*12)
end

%Filter to treat aliasing. Evaluate each pixel separately, if the same
%pixel in the following frame makes a sudden "jump" in the value, the
%pixels is aliased and are set to the real value before next time step is
%evaluated. The limitation is that if the first frame is aliased,
%everything will be fucked up. To solve this, the first frame is filtered
%by comparing the suspicious pixels with the mean velocity of all pixels
%in the ROI, and apply a threshold to decide wether the pixel is aliasied
%or not. The following frames are filtered by comparing the pixels with the
%former frame.
%OBS! Subject 22 is seriously aliased with values that exceeds 2*V_enc,
%therefore a third filter is applied which sets these values to +-2*V_enc.

teller = 1;

th_min = -3.5;%-3.5;
th_max = 3.5;%3.5;
th_al = 5.5;
%OBS! Special treatment for Subject_22 due to psyko aliasing. First frame
%need to be filtered like: for i=126:svv(3) and j=155:svv(3) with
%th_min/max -2/2 and th_al 5.2
%Also Subject 29, best result with NO filter.

for k = 1
    for i=1:svv(2)
        for j = 1:svv(2);
            if abs(peakvel(k,1))> th_al && abs(TMP(i,j,k)) ~=0 || abs(peakvel(k,2))> th_al && abs(TMP(i,j,k)) ~=0 
                if vk(k)>0 && TMP(i,j,k)<th_min
                    TMP(i,j,k) = TMP(i,j,k) + 2*V_enc; 
                end
                if vk(k)<0 && TMP(i,j,k)>th_max
                    TMP(i,j,k) = TMP(i,j,k) - 2*V_enc; 
                 teller = teller+1;   
                end
            end
            
        end
    end
end


for k = 2:svv(3);
    for i = 1:svv(2)
        for j = 1:svv(2)
            if abs(TMP(i,j,k) - TMP(i,j,k-1)) > 1.1*V_enc && abs(TMP(i,j,k)) ~=0
                if abs(TMP(i,j,k))/TMP(i,j,k) > 0
                    TMP(i,j,k) = TMP(i,j,k) - 2*V_enc;
                elseif TMP(i,j,k) ~=0
                   TMP(i,j,k) = TMP(i,j,k) + 2*V_enc;
                end
            end
        end
    end
end
% % 
% %%Additional filter for subject_22
% for k = 2:svv(3);
%     for i = 1:svv(2)
%         for j = 1:svv(2)
%             if abs(TMP(i,j,k) - TMP(i,j,k-1)) > 1.9*V_enc && abs(TMP(i,j,k)) ~=0
%                 if abs(TMP(i,j,k))/TMP(i,j,k) > 0
%                     TMP(i,j,k) = -2*V_enc;
%                 elseif TMP(i,j,k) ~=0
%                    TMP(i,j,k) = 2*V_enc;
%                 end
%             end
%         end
%     end
%     maxmin(k,1) = max(max(TMP(:,:,k).*double(mask_image)));
%     maxmin(k,2) = min(min(TMP(:,:,k).*double(mask_image)));
%end

% 
for k = 1:svv(3);%points%start:start+5;
    figure;imagesc(TMP(:,:,k));colorbar;colormap('jet');
    title(['Filter', num2str(k)]);
    axis(xylim);caxis([-1 1]*12)
end

%Next: Extract mean velocity (vk) and volumetric flow rate (Q), where 
%flow rate includes total, positive part and negative part.
sumP = zeros(svv(3),1);
sumN = sumP;
tellerp = sumP;
tellern = sumP;

for k = 1:svv(3)
    vk(k) = sum(sum(TMP(:,:,k).*mask_image));
    tmp2 = TMP(:,:,k);
    Q(k,1) = sum(tmp2(idx))*dx*dy;
    for j = 1:length(idx)
        if tmp2(idx(j))>0
            sumP(k) = sumP(k) + tmp2(idx(j));
            tellerp(k) = tellerp(k)+1;
        elseif tmp2(idx(j))<0
            sumN(k) = sumN(k) + tmp2(idx(j));
            tellern(k) = tellern(k)+1;
        end
    end
     Q(k,2) = sumP(k)*dx*dy;
     Q(k,3) = sumN(k)*dx*dy;
end

vk = vk/numpix;

%Qint is volumetric flow rate integrated over dt, which is the way Alperine
%define V.

Qint = Q(:,1)*dt;

%Pressure gradient 
dp = zeros(svv(1), svv(2), svv(3));
dpa = zeros(svv(1), svv(2), svv(3));
dpmu = zeros(svv(1), svv(2), svv(3));

for k=1:svv(3)-1 
  for i = 2:svv(1)-1;
    for j = 2:svv(1)-1;  
      if mask_image(i,j) > 0 
       dp(i,j,k) = (-density*(1.0*TMP(i,j,k+1) - TMP(i,j,k)) / dt) +...
           visc*((TMP(i+1,j,k)-2*TMP(i,j,k)+TMP(i-1,j,k)/(dx.^2) + (TMP(i,j+1,k)-2*TMP(i,j,k)+TMP(i,j-1,k))/(dy.^2))); 
       dpa(i,j,k) = (-density*(1.0*TMP(i,j,k+1) - TMP(i,j,k)) / dt);
       dpmu(i,j,k) =  visc*((TMP(i+1,j,k)-2*TMP(i,j,k)+TMP(i-1,j,k)/(dx.^2) + (TMP(i,j+1,k)-2*TMP(i,j,k)+TMP(i,j-1,k))/(dy.^2))); 
      end
    end
  end
end

%Unit of the calculated dp/dz = (1/10*Pa) /cm
%Wish mmHg, convert: 1 mmHg = 133.322 Pa -> 
%ans (mmHg) = ans (0.1 Pa)*(1/133.322 mmHg/Pa).
dp_mean = sum(sum(dp)); 
dp_mean = reshape(dp_mean, [1 length(dp_mean)]);
dp_mean = dp_mean / numpix;
dp_mean = dp_mean / 1333.2;

%Save to file

% fid = fopen(['20151022',subject_folders(p).name,'_VolumetricCSFFlow'],'w');
% for i = 1:svv(3)
%     fprintf(fid, '%f %f %f %f\n',time(i),Q(i,1),Q(i,2),Q(i,3));
% end
% fclose(fid);
% 
% % fid2 = fopen([subject_folders(p).name,'_VolumetricCSFFlowInt'],'w');
% % for i = 1:svv(3)
% %     fprintf(fid2, '%f %f\n',time(i),Qint(i));
% % end
% % fclose(fid2);
% % 
% fid3 = fopen(['20151022',subject_folders(p).name,'_MeanVelCSF'],'w');
% for i = 1:svv(3)
%     fprintf(fid3, '%f %f\n',time(i),vk(i));
% end
% fclose(fid3);
% 
% fid4 = fopen(['20151022',subject_folders(p).name,'_PressGradCSF'],'w');
% for i = 1:svv(3)
%     fprintf(fid4, '%f %f\n',time(i),dp_mean(i));
% end
% fclose(fid4);
% 
% fid5 = fopen(['20151022',subject_folders(p).name,'_nbPosNeg'],'w');
% for i = 1:svv(3)
%     fprintf(fid5, '%f %f %f %f\n',time(i),numpix,tellerp(i),tellern(i));
% end
% fclose(fid5);

figure
plot(time,vk,'ko-')
legend('Meanvel')
xlabel('Time  [s]')
ylabel('Mean velocity  [cm/s]')
title(['Subject\_',subject_folders(p).name,'   CSF'])


figure
plot(time,Q(:,1),'ko-',time,Q(:,2),'r*-',time,Q(:,3),'b*-')
xlabel('Time   [s]')
ylabel('Volumetric flow rate   [ml/s]')
title(['Subject\_',subject_folders(p).name,'   CSF'])

figure
plot(time,Qint,'mp-')
xlabel('Time   [s]')
ylabel('Volumetric flow    [ml]')
title(['Subject\_',subject_folders(p).name,'   CSF'])

figure
plot(time(1:end-1), dp_mean(1:end-1),'ko-','LineWidth',2,'MarkerFaceColor','k')
axis([0 1.1 -0.015 0.015])
xlabel('Time  [s]')
ylabel('dP  [mmHg/cm]')
%title(['Subject\_',subject_folders(p).name,'   CSF'])

figure
plot(time,(tellerp/numpix)*100,'r-o',time,(tellern/numpix)*100,'b-o')
xlabel('Time   [s]')
ylabel('% pos/neg')
title(['Subject\_',subject_folders(p).name,'   CSF'])

% ci = 23;
% figure;imagesc(velocity_image(:,:,ci))
% colormap(gray)
% colorbar
% hold on
% plot(xx(idx),yy(idx),'r.');
% title(subject_folders(p).name)


end