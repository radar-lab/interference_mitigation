close all;
clear;
clc;

%% Provide the filename
% Radar raw data captured by TI AWR1243BOOST + TI mmwave DevPack + TI TSW1400
% TI AWR1243BOOST:      http://www.ti.com/tool/AWR1243BOOST
% TI mmwave DevPack:    http://www.ti.com/tool/MMWAVE-DEVPACK  
% TI TSW1400:           http://www.ti.com/tool/TSW1400EVM
% Configuration:        https://training.ti.com/sites/default/files/docs/mmWave_sensor_raw_data_capture_using_TSW1400_board_v3.pdf
fname = '20181202-193712.bin'; % Fig. 6: Increased SIR
%fname = '20181128-214451.bin'; % Fig. 7: Ghost Target 

%% Waveform Parameters
Fc          = 77e9;         % Carrier frequency
c           = 3e8;          % Speed of light
BW          = 750e6;        % Chirp bandiwdth
PRI         = 29.56e-6;     % Chirp duration
ChirpRate   = 29.306e12;    % Chirp rate
ADCFs       = 20e6;         % ADC sampling rate
NFast       = 512;          % Number of ADC samples during one chirp
NPRI        = 128;          % Number of chirps during one CPI
n_adc_bits  = 16;           % NUmber of the ADC bits
n_rx        = 4;            % Number of receiving antenna channel
n_tx        = 1;            % Number of transmitting antenna channel

%% Open file and pre-porcess the data
fid     = fopen(fname,'r');     % Open the file in read mode
dat     = fread(fid, NFast * NPRI * n_rx * n_tx * 2,'uint16'); % Read the data into Matlab enviroment
fclose(fid);                    % Close the file
% Form the data
dat     = dat - 2^15;
dat     = reshape(dat,8,[]); 
dat     = dat([1,3,5,7,2,4,6,8],:);
cdat    = dat(1:4,:) + 1i*dat(5:8,:);
norm_factor = (sqrt(2)/(2^(n_adc_bits-1)-1));
cdat    = cdat*norm_factor;
% Form the data for each antenna channel
RX1     = cdat(1,:);
RX2     = cdat(2,:);
RX3     = cdat(3,:);
RX4     = cdat(4,:);
% Reshape the data
RX1     = reshape(RX1, NFast, NPRI);
RX2     = reshape(RX2, NFast, NPRI);
RX3     = reshape(RX3, NFast, NPRI);
RX4     = reshape(RX4, NFast, NPRI);
% Only use the data from receiver #1 and disregard the other channel's data
adc_dat = RX1; 

figure
%% Test
RangeDoppler = (adc_dat);
RangeDoppler = fft(RangeDoppler, NFast, 1); % FFT for each colum
RangeDoppler = (abs(fftshift(fft(RangeDoppler, NPRI, 2)))); % FFT for each row
dBRangeDopplerScaled = db(abs(RangeDoppler));
dBRangeDopplerScaled(dBRangeDopplerScaled<20) = 20;
[Doppler, Range] = meshgrid(-NPRI/2:1:NPRI/2-1, -NFast/2:1:NFast/2-1);
% subplot(122);
surf(Doppler, Range, dBRangeDopplerScaled, 'LineStyle', 'none');
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('(e)');
xlabel('Doppler FFT bin');
xlim([-NPRI/2 NPRI/2-1]);
ylabel('Range FFT bin');
ylim([-NFast/2 NFast/2-1]);
zlabel('Power(dB)');
% zlim([20 55]);
colorbar;
view([58,40]); 

figure
%% Apply adaptive noise canceler on one chirp data
OneChirp                = adc_dat(:, 100);  % Select the 50th chirp
OneRangeFFT             = fftshift(fft(OneChirp, NFast));
set(gcf, 'Position', [0 0 1000 800])
subplot(411)
plot(real(OneChirp));
xlim([1 NFast]);
ylim([-0.7 0.7]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
title('(a)');
xlabel(['Time (' 956 's)']);
ylabel('Amplitude(linear)');
xticklabels({'2.5', '5', '7.5', '10', '12.5', '15', '17.5', '20', '22.5', '25'});
subplot(412)
FFTBin = (-NFast/2:1:NFast/2-1);
plot(FFTBin, db(abs(OneRangeFFT)));
xlim([-NFast/2 NFast/2-1]);
ylim([-20 60]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
title('(b)');
xlabel('FFT Bin');
ylabel('Amplitude(dB)');
subplot(413)
FFTBin = (-NFast/2:1:NFast/2-1);
plot(FFTBin, (real(OneRangeFFT)));
xlim([-NFast/2 NFast/2-1]);
% ylim([-20 20]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
title('(c)');
xlabel('FFT Bin');
ylabel('Amplitude(linear)');
subplot(414);
FFTBin = (-NFast/2:1:NFast/2-1);
plot(FFTBin, (imag(OneRangeFFT)));
xlim([-NFast/2 NFast/2-1]);
% ylim([-20 20]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
title('(d)');
xlabel('FFT Bin');
ylabel('Amplitude(linear)');