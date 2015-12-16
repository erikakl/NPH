close all

files = dir('Erika*_VolumetricCSFFlow');
numfiles = length(files);
Q = zeros(numfiles,3);

for k = 1:numfiles
    a = load(files(k).name);
    figure
    plot(a(:,1),a(:,2),'k*-',a(:,1),a(:,3),'r*-',a(:,1),a(:,4),'b*-')
    title(['Erika ',num2str(k)])
    legend('Q-tot','Q-pos','Q-neg')
    xlabel('Time [s]')
    ylabel('Volumetric flow rate [ml/s]')
    Q(k,1) = trapz(a(:,1),a(:,2));
    Q(k,2) = trapz(a(:,1),a(:,3));
    Q(k,3) = trapz(a(:,1),a(:,4));
    
end

fid = fopen('Erika_Flux','w');
for i = 1:numfiles
    fprintf(fid,'%f %f %f\n',Q(i,1),Q(i,2),Q(i,3));
end

b = load('Erika_Flux');
disp(b)


