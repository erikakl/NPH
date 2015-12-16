a1 = load('Erika02_Peakvel_Av');
b1 = load('Erika05_Peakvel_Av');

a2 = load('Erika02_MeanVelCSF');
b2 = load('Erika05_MeanVelCSF');

a3 = load('Erika02_VolumetricCSFFlow');
b3 = load('Erika05_VolumetricCSFFlow');

a4 = load('Erika02_PressGradCSF');
b4 = load('Erika05_PressGradCSF');


figure
plot(a1(:,1),a1(:,2),b1(:,1),b1(:,2))
title('Peak velocity')
xlabel('Time [s]')
ylabel('Velocity [cm/s]')
legend('T1 PPU2','T2 PPU2')


