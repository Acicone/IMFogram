%% Example of Twinkle Twinkle Little Star Piano song - first four notes Do - Do - Sol - Sol
%  Please cite:
%
%  P. Barbe, A. Cicone, W. S. Li, H. Zhou. 
%  "Time-frequency representation of nonstationary signals: the IMFogram". 
%  Pure and Applied Functional Analysis, 2021. 
%  arXiv http://arxiv.org/abs/2011.14209

clear
close all
clc

[y,fs] = audioread('Little_star_by_Alice.mp4');

figure
plot(y(fs*21:fs*23))

sound(y(fs*21:fs*23),fs)
y=y(fs*21:fs*23);

%% Spectrogram
figure
tic
pspectrum(y,fs,'spectrogram','FrequencyLimits',[0,2000])
time_Spectrogram=toc;
title('')
hold on
plot3(0.7959,0.523,1,'xr','markersize',40,'linewidth',5)
plot3(0.22,0.523,1,'xr','markersize',40,'linewidth',5)
plot3(1.265,0.784,1,'xk','markersize',40,'linewidth',5)
plot3(1.82,0.784,1,'xk','markersize',40,'linewidth',5)
legend('Do = 523 Hz','Do = 523 Hz','Sol = 784 Hz','Sol = 784 Hz')
axis tight
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'fontsize', 30);

%% Spectrogram in 3D
figure
pspectrum(y,fs,'spectrogram','FrequencyLimits',[0,2000])
title('')
axis tight
view(330,67)
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'fontsize', 30);

%% FIF decomposition into IMFs

opts=Settings_FIF_v3('alpha',30,'NumSteps',5,'ExtPoints',20,'NIMFs',15,'verbose',0,'alpha','Median');
[IMFs,stats]=FIF_v2_12(y,opts);

%% IMFs plotting

plot_imf_v10([y';IMFs],(0:fs*2)/fs,9)
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'fontsize', 30);

%% IMFogram

MaxF=2000;
winLen=2000;
winOverPerc=90;
NiF=128;

IMFogram_v1(IMFs(:,1:end),fs,winOverPerc,[0 MaxF],NiF,winLen);
toc
hold on
plot3(0.7959,523,1,'xr','markersize',40,'linewidth',5)
plot3(0.22,523,1,'xr','markersize',40,'linewidth',5)
plot3(1.265,784,1,'xk','markersize',40,'linewidth',5)
plot3(1.82,784,1,'xk','markersize',40,'linewidth',5)
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'fontsize', 30);


%% IMFogram in 3D

IMFogram_v1(IMFs(:,1:end),fs,winOverPerc,[0 MaxF],NiF,winLen);
view(330,67)
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'fontsize', 30);

