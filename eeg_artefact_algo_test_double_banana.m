%.........................................................................................%
% 2021-08-24 Created by Pathan Fayaz Khan, Signal Processing lead, Neurostellar.
% This program removes eyeblink, Eye movement, Electrode Pop, Baseline drift and 
% Muscle artefact corrections using a hybrid combination of Canonical
% Correlation Analysis and Ensemble Emperical mode Decomposition
%..........................................................................................%

%.....................FILE HANDLING.............................%
clc; clear all; close all;
ft_defaults
nicoletfile = 'Patient234_DESKTOP-CRU7CPI_t1.e';
hdr = ft_read_header(nicoletfile);
dat = ft_read_data(nicoletfile, 'header', hdr);
dat = dat';
%.....................END OF FILE HANDLING......................%

%.....................CREDENTIALS...............................%
%% 
fs = hdr.Fs;                         %Sampling Frequency
seg_length = 10;                     %Length of each segment to be processed
ResRequired = 0.1;
nfft= fs/ResRequired;                                                       % number of overlapping time samples
window= 256;
noverlap= window/2;
[nr, nc] = size(dat);
for i = 1:nc
    dat_cor(:,i) = dat(seg_length*fs:nr,i) - mean(dat(seg_length*fs:nr,i));    %Mean correction of Raw data
end
nr5 = length(dat_cor);
channel_labels = {'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', 'Fz', 'P3', 'P4', 'Pz', 'C3', 'C4', 'Cz', 'O1', 'O2', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'ECG', 'Photic'};
[nl, cl] = size(channel_labels);
q = 1; dat_corr = [];
for i = 1:cl
    [logic(i, 1), channel_numb(i,1)] = ismember(channel_labels(1,i),hdr.label);
    if(logic(i,1)==1)
        dat_corr(1:nr5,q) = dat_cor(:,channel_numb(i,1));
        channel_final_labels(q,1) = channel_labels(1,i); 
        q = q+1;
    end
end
%.....................END OF CREDENTIALS........................%
%.....................NOTCH FILTERS.............................%
%% 
w1 = 50/(fs/2); 
w2 = 100/(fs/2);
bw1 = w1/35;
bw2 = w2/35;
[b1,a1] = iirnotch(w1,bw1);
[b2,a2] = iirnotch(w2,bw2);
eeg_n50 = filtfilt(b1,a1,dat_corr);
eeg_n100 = filtfilt(b2,a2,eeg_n50);
%.....................END OF NOTCH FILTERS......................%

%.....................FILTERING.................................%
f1 = 0.3; f2 = 70;
wn1 = f1/(fs/2); wn2 = f2/(fs/2);
wn = [wn1,wn2];
[B,A] = butter(2, wn);
eeg_bp = filtfilt(B,A,eeg_n100);
%.....................END OF FILTERING..........................%

%.....................MAIN PROGRAM..............................%                                                                    %No. channels
%% 
[nr6, nc6] = size(eeg_bp);
M = seg_length*fs;                                                                  %Length of each segment
N3 = round(nr/M);                                                           %No. of 20 sec segments
cleaned_sig = [];
NumSeg = 150;
eeg_Raw = eeg_n100(1:NumSeg*seg_length*fs,:);
ecg_photic = 0;
if ((logic(length(logic)-1,1)==1))
    ecg_photic = ecg_photic+1;
end
if ((logic(length(logic),1)==1))
    ecg_photic = ecg_photic+1;
end
for i = 1:NumSeg
    eeg_seg(1:M,:) = eeg_bp((M*(i-1))+1:(M*i),1:nc6-ecg_photic);
    eeg_seg_Delay = delayseq(eeg_seg,3);
    [A, B, r, U, V] = canoncorr(eeg_seg, eeg_seg_Delay);
    U4 = U;
    [nr1, nc1] = size(U);
    for p = 1:nc1
        if(r(1,p)<0.4)
            U(:,p)=zeros(nr1,1);
        end
    end
    U1 = U';
    for j = 1:nc1
        [imf] = emd_sssc(U1(j,:),fs, 'display', 0);
        imf1 = imf';
        [nr2, nc2] = size(imf);
        imf_clean = zeros(M,1);
        for k = 1:nr2
            [pxx1,~]= pwelch(imf1(:,k),window,noverlap,nfft,fs,'onesided','power');
            coefs1_power = mean(pxx1(1:20,1)) * 20;
            coefs2_power = mean(pxx1(21:300,1)) * 300;
            ratio = coefs2_power/coefs1_power;
            disp('Analyzing....')
            if(ratio>1)
              imf_clean = imf_clean+imf1(:,k);
            end
        end
        U(:,j) = imf_clean;
    end
    cleaned_sig_temp = U*pinv(A);
    cleaned_sig = vertcat(cleaned_sig, cleaned_sig_temp);
end
cleaned_sig_final = horzcat(cleaned_sig,dat_corr(1:NumSeg*seg_length*fs,nc6-1:nc6));
%..................END OF MAIN PROGRAM....................%

%% 
%......................PLOTTING...........................%

% channel = 1;
% t = 0:1/fs:NumSeg*seg_length-(1/fs);
% close all;
% figure;
% subplot(2,1,1)
% plot(t, eeg_Raw(:,channel)); 
% title('CCA EEMD HYBRID ALGORITHM')
% xlabel('Time (s)')
% ylabel('Amplitude (microVolt)')
% legend('Raw EEG');
% % xlim([0 20])
% ylim([-200 200])
% grid on;
% 
% subplot(2,1,2)
% plot(t, cleaned_sig_final(:,channel)); 
% xlabel('Time (s)')
% ylabel('Amplitude (microVolt)')
% legend('Cleaned EEG');
% % xlim([0 20])
% ylim([-200 200])
% grid on;

%.............POWER SPECTRAL DENSITY...............% 
% [pxx1,fp1]= pwelch(cleaned_sig_final(:,channel),window,noverlap,nfft,fs,'onesided','power');
% [pxx2,fp2]= pwelch(eeg_Raw(:,channel),window,noverlap,nfft,fs,'onesided','power'); 
% figure;
% plot(fp1,10*log10(pxx1)); 
% hold on; 
% plot(fp2,10*log10(pxx2),'r');
% title('CCA EEMD HYBRID ALGORITHM')
% xlabel('Frequency (Hz)')
% ylabel('Power (microVolt^2)')
% legend('Cleaned EEG','Raw EEG');
% % xlim([0 20])
% % ylim([-100 100])
% grid on;
% %...................FINAL END......................%
% % plot(t, eeg_seg(:,channel)-cleaned_sig(:,channel));

save cleaned_sig_final_234.txt cleaned_sig_final '-ascii';
save eeg_Raw_234.txt eeg_Raw '-ascii';
writecell(channel_final_labels,'channel_final_labels_234.txt');
save fs.txt fs '-ascii';

% cleaned_sig_final = textread('cleaned_sig_final_715.txt');
% eeg_Raw = textread('eeg_Raw_715.txt');
% channel_final_labels = readcell('channel_final_labels_715.txt');
% fs = textread('fs.txt');
% 
% [nr8, nc8] = size(cleaned_sig_final);
% mixed_sig = horzcat(eeg_Raw(:,1:round(nc8/2)), cleaned_sig_final(:, 1:round(nc8/2)));
% figure;
% [nr, nc] = size(mixed_sig);
% k = 1;
% amp = 1;
% fig = gcf;
% timeBase = 15*fs; % 15 sec
% XTicks = 0:fs:15*fs; %15 sec
% YTicks = -200:200:(nc)*200;
% [nr9, nc9] = size(channel_final_labels);
% channel_final_labels_half_numb = channel_final_labels(1:round(nr9/2),1);
% YLabel = [{''},channel_final_labels_half_numb',channel_final_labels_half_numb', {''}];
% YLabel = flip(YLabel');
% while(k<=nr)
%   for i = 1:nc
%      plot((200*(i-1)) + (amp*mixed_sig(k:k+timeBase,(nc+1)-i)));
%      ax = fig.CurrentAxes;
%      ax.XLim = [0 timeBase];
%      ax.YLim = [-200 nc*200];
%      ax.XTick = XTicks;
%      ax.YTick = YTicks;
%      ax.XTickLabel = round((k/fs):1:((k+timeBase)/fs)); 
%      ax.YTickLabel = YLabel;
%      ax.FontSize = 14; ax.FontWeight = 'bold'; ax.XGrid = 'on'; ax.YGrid = 'on';
%      hold on;
%   end
%   was_a_key = waitforbuttonpress;
%   if ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'rightarrow')))
%       k = k+500;
%     if(k>(length(eeg_Raw)-timeBase))
%         k = length(eeg_Raw)-timeBase;
%     end
%   elseif ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'leftarrow')))
%     k = k-500;
%     if(k<1)
%         k = 1;
%     end
%   elseif ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'uparrow')))
%       amp = amp+0.1;
%       if(amp>5)
%           amp = 5;
%       end
%   elseif ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'downarrow')))
%       amp = amp-0.1;
%       if(amp<1)
%           amp = 1;
%       end
%   end
%    hold off; 
% end
% 
% %...................................................................%
% banana_ch = {'Fp1-F3', 'F3-C3', 'C3-P3', 'P3-O1', 'Fp2-F4', 'F4-C4', 'C4-P4', 'P4-O2', 'Fp1-F7', 'F7-T1', 'T1-T3', 'T3-T5', 'T5-O1', 'Fp2-T8', 'F8-T2', 'T2-T4', 'T4-T6', 'T6-O2', 'Fz-Cz', 'Cz-Pz', 'ECG', 'Photic'};
% [nr6, nc6] = size(channel_final_labels);
% [nr7, nc7] = size(cleaned_sig_final);
% w = 1;
% ban_config_ch = [];
% 
% [v1, u1] = ismember({'Fp1','F3'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%    ban_config_label(w,1) =  banana_ch(1,1);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'F3','C3'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%    ban_config_label(w,1) =  banana_ch(1,2);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'C3','P3'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%    ban_config_label(w,1) =  banana_ch(1,3);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'P3','O1'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%    ban_config_label(w,1) =  banana_ch(1,4);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'Fp2','F4'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%    ban_config_label(w,1) =  banana_ch(1,5);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'F4','C4'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,6);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'C4','P4'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,7);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'P4','O2'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,8);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'Fp1','F7'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,9);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'F7','T1'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,10);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'T1','T3'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,11);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'T3','T5'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,12);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'T5','O1'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,13);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'Fp2','T8'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,14);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'F8','T2'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,15);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'T2','T4'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,16);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'T4','T6'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,17);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'T6','O2'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,18);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'Fz','Cz'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,19);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'Cz','Pz'},channel_final_labels);
% if(sum(v1)==2)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1))-cleaned_sig_final(:,u1(1,2));
%       ban_config_label(w,1) =  banana_ch(1,20);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'ECG'},channel_final_labels);
% if(sum(v1)==1)
% %    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1));
%    ban_config_ch(1:nr7,w) = zeros(nr7,1);
%       ban_config_label(w,1) =  banana_ch(1,21);
%    w = w+1;
% end
% 
% [v1, u1] = ismember({'Photic'},channel_final_labels);
% if(sum(v1)==1)
%    ban_config_ch(1:nr7,w) = cleaned_sig_final(:,u1(1,1));
%       ban_config_label(w,1) =  banana_ch(1,22);
%    w = w+1;
% end
% 
% [nr, nc] = size(ban_config_ch);
% k = 1;
% amp = 1;
% fig = gcf;
% timeBase = 15*fs; % 15 sec
% XTicks = 0:fs:15*fs; %15 sec
% YTicks = -200:200:(nc)*200;
% YLabel = [{''},ban_config_label', {''}];
% YLabel = flip(YLabel');
% 
% while(k<=nr)
%   for i = 1:nc
%      plot((200*(i-1)) + (amp*ban_config_ch(k:k+timeBase,(nc+1)-i)));
%      ax = fig.CurrentAxes;
%      ax.XLim = [0 timeBase];
%      ax.YLim = [-200 nc*200];
%      ax.XTick = XTicks;
%      ax.YTick = YTicks;
%      ax.XTickLabel = round((k/fs):1:((k+timeBase)/fs)); 
%      ax.YTickLabel = YLabel;
%      ax.FontSize = 14; ax.FontWeight = 'bold'; ax.XGrid = 'on'; ax.YGrid = 'on';
%      hold on;
%   end
%   was_a_key = waitforbuttonpress;
%   if ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'rightarrow')))
%       k = k+500;
%     if(k>(length(eeg_Raw)-timeBase))
%         k = length(eeg_Raw)-timeBase;
%     end
%   elseif ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'leftarrow')))
%     k = k-500;
%     if(k<1)
%         k = 1;
%     end
%   elseif ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'uparrow')))
%       amp = amp+0.1;
%       if(amp>5)
%           amp = 5;
%       end
%   elseif ((was_a_key)&&(strcmp(get(fig,'CurrentKey'), 'downarrow')))
%       amp = amp-0.1;
%       if(amp<1)
%           amp = 1;
%       end
%   end
%    hold off; 
% end
