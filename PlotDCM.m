close all
clear all

%for nm = 10

%info = dicominfo([num2str(nm),'.dcm']);
info = dicominfo('51#001.dcm');
Y = dicomread(info);
figure, imshow(Y);
imcontrast;

figure
imagesc(Y)
colorbar
title('Subject 18');

disp(max(max(Y)))
disp(min(min(Y)))

%end