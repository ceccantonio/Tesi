clear all
close all

buff = readmatrix('gamma_sweep.txt');

gamma = unique(buff(:,1));
f = unique(buff(:,2));
d = reshape(buff(:,3),length(f),length(buff(:,3))/length(f));

%%
close all

figure(1)
set(gcf,'Position',[650 250 300 600])
hold on
c = winter(size(d,2));
for j = 1:size(d,2)
    plot(f*1e-9,log(d(:,j))+3*j,'Linewidth',3,'Color',c(j,:))
end
grid on
axis tight
xlabel('Frequency [GHz]')

%%
figure(2)
pcolor(f*1e-9,gamma,log(d'))
shading flat
colormap turbo
xlabel('Frequency [GHz]')
ylabel('Angle [rad]')
