close all
clear all


addpath C:\Users\erikakl\Documents\MATLAB\spm12\

folder_subjects = strcat(pwd,'/ErikaTime2');
subject_folders = dir(folder_subjects);
subject_folders = subject_folders(3:end);

velocity_name = 'ccj.nii';
mask_name = 'ccj_mask.nii';
dicom_name = 'hdr.dcm';

%Loop for all subjects.
for p = 2
  
    
current_folder = strcat(folder_subjects,'/',subject_folders(p).name);
file_status_ccj = dir(fullfile(current_folder,velocity_name));                           
velocity_volume = spm_vol(fullfile(current_folder,velocity_name));
mask_volume = spm_vol(fullfile(current_folder,mask_name));

density = 1.0007; % density for CSF is 1.0007 g/(cm^3), (this value comes from Alperin)
visc = 1.1*0.01;%mu for CSF is about 1.1 cP ("centipoise"), where 1 cP = 0.01 g/(cm*s) (this value comes from Alperin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find nb mask
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for N = 1:velocity_volume.dim(3)
    mask_image = flipdim(spm_slice_vol(mask_volume,spm_matrix([0 0 N]),mask_volume.dim(1:2),1),2);
    figure;imagesc(mask_image);
    title(num2str(N))
    axis([120 160 130 190]);
end
plotend = velocity_volume.dim(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frame no Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p == 1 
N = 18;
else N = 17;
end

velocity_image = zeros(velocity_volume.dim(1),velocity_volume.dim(2),velocity_volume.dim(3)); 
mask_image = flipdim(spm_slice_vol(mask_volume,spm_matrix([0 0 N]),mask_volume.dim(1:2),1),2);
figure;imagesc(mask_image);
title('mask image Org')
axis([100 200 100 200]);
%info = dicominfo(fullfile(current_folder,dicom_name));
%V_enc = info.Private_2001_101a(3);%Max value possible, if exceeded: aliasing
svv = size(velocity_image);
tmpMask = mask_image;
TestMask = tmpMask;
for i = 1:svv(1)
    for j = 2:svv(2)
        if tmpMask(i,j) == 0;
           tmpMask(i,j) = 0;
        else tmpMask(i,j) = 1;
        end
        
    end
end


tmpMask(1:288,1)=0;
tmpMask(1:288,288)=0;
tmpMask(1,1:288)=0;
tmpMask(288,1:288)=0;

%tmpmask = double(mask_image);
%tmpmask = ((2*V_enc)/4.095)*mask_image(:,:)-V_enc;
%TMPmask = tmpmask;%(tmp.*double(mask_image));     
figure
imshow(tmpMask)
title('background 0')
%Convert ROI to black and white.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%REMEMBER!!!
mask_image = im2bw(tmpMask,1e-6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NB!!!!! Fliplr 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask_image = fliplr(mask_image);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove noisy pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p == 3
    mask_image(140,146,:) = 0;
    mask_image(171,143,:) = 0;
end

%figure;plot(mask_image,'.');title('Flipped')
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
disp(max(max(max(velocity_image))))
disp(min(min(min(velocity_image))))

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


%convert pix-values to cm/s
%6=6cm/s (V_enc), 2.0475=0cm/s, -6=-6cm/s(-V_enc)
for k=1:svv(3)
    tmp = double(velocity_image(:,:,k));
  %  tmp = ((2*V_enc)/4.095)*tmp-V_enc;
    vk(k) = sum(tmp(idx));
    peakvel(k,1) = max(max(tmp(idx)));
    peakvel(k,2) = min(min(tmp(idx)));
end


%Compute mean velocity, average over number of pixels
vk = vk/(length(idx));


points =1:svv(3);

%For animations, zoom ROI for each subject
%aviobj = VideoWriter('mymovie.avi'); 
xylim = [100  200 100 200];


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
%[130 180 120 170] %51



TMP = zeros(svv(1),svv(2),svv(3));
for k=1:svv(3)
    tmp = double(velocity_image(:,:,k));
%    tmp = ((2*V_enc)/4.095)*tmp-V_enc;
    TMP(:,:,k) = (tmp.*double(mask_image));     
end


maxmin = zeros(svv(3),2);

%Plot for safety check.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

for k = 1%:svv(3);%points
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

th_min = -3.5;
th_max = 3.5;
th_al = 5.5;
%OBS! Special treatment for Kent03 due to psyko aliasing. First frame
%need to be filtered like: with
%th_min/max -2/2 and th_al 5.2
%OBS! In general th_al = 5.5 (for all subjects except 22 and 28)
%28: th_al = 5.2


for k = 1
    for i=1:svv(2)
        for j = 1:svv(2);
            if abs(peakvel(k,1))> th_al && abs(TMP(i,j,k)) ~=0 || abs(peakvel(k,2))> th_al && abs(TMP(i,j,k)) ~=0 
                if vk(k)>0 && TMP(i,j,k)<th_min
                    TMP(i,j,k) = TMP(i,j,k) + 2*V_enc; 
                end
                if vk(k)<0 && TMP(i,j,k)>th_max
                    TMP(i,j,k) = TMP(i,j,k) - 2*V_enc; 
                    
                end
            end
            
        end
    end
    maxmin(k,1) = max(max(TMP(:,:,k)));
    maxmin(k,2) = min(min(TMP(:,:,k)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For p=2,3,5 use constant 1.9 instead of 1.1
% Earilier: if abs(TMP(i,j,k) - TMP(i,j,k-1)) > 1.1*V_enc && abs(TMP(i,j,k)) ~=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
for k = 2:svv(3);
    for i = 1:svv(2)
        for j = 1:svv(2)
            if abs(TMP(i,j,k) - TMP(i,j,k-1)) > 1.9*V_enc && abs(TMP(i,j,k)) ~=0
                if abs(TMP(i,j,k))/TMP(i,j,k) > 0
                    TMP(i,j,k) = TMP(i,j,k) - 2*V_enc;
                elseif TMP(i,j,k) ~=0
                   TMP(i,j,k) = TMP(i,j,k) + 2*V_enc;
                end
            end
        end
    end
    maxmin(k,1) = max(max(TMP(:,:,k).*double(mask_image)));
    maxmin(k,2) = min(min(TMP(:,:,k).*double(mask_image)));
end



for k = 1%:svv(3);%plotend; %start:start+5;
    figure;imagesc(TMP(:,:,k));colorbar;colormap('jet');
    title(['Filter', num2str(k)]);
    axis(xylim);caxis([-1 1]*12)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find min/max 10% for each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb = round(numpix*0.1);
storepeak = zeros(svv(3),2);
for k = 1:svv(3)
    velvec = zeros(numpix,1);
    tmpv = TMP(:,:,k);
    for i = 1:numpix
        velvec(i) = tmpv(idx(i));
      
    end
    velvec2 = sort(velvec);
    storepeak(k,1) = mean(velvec2(end-nb:end));
    storepeak(k,2) = mean(velvec2(1:nb));
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

% fid = fopen(['Erika',subject_folders(p).name,'_VolumetricCSFFlow'],'w');
% for i = 1:svv(3)
%     fprintf(fid, '%f %f %f %f\n',time(i),Q(i,1),Q(i,2),Q(i,3));
% end
% fclose(fid);
% 
% fid2 = fopen(['Ny_Erika',subject_folders(p).name,'_Peakvel_Av'],'w');
% for i = 1:svv(3)
%     fprintf(fid2, '%f %f %f\n',time(i),storepeak(i,1),storepeak(i,2));
% end
% fclose(fid2);

% 
% fid3 = fopen(['Erika',subject_folders(p).name,'_MeanVelCSF'],'w');
% for i = 1:svv(3)
%     fprintf(fid3, '%f %f\n',time(i),vk(i));
% end
% fclose(fid3);
% 
% fid4 = fopen(['Erika',subject_folders(p).name,'_PressGradCSF'],'w');
% for i = 1:svv(3)
%     fprintf(fid4, '%f %f\n',time(i),dp_mean(i));
% end
% fclose(fid4);

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
plot(time, dp_mean,'bp-')
xlabel('Time  [s]')
ylabel('dP  [mmHg/cm]')
title(['Subject\_',subject_folders(p).name,'   CSF'])

figure
plot(time,(tellerp/numpix)*100,'r-o',time,(tellern/numpix)*100,'b-o')
xlabel('Time   [s]')
ylabel('% pos/neg')
title(['Erika   CSF'])

end