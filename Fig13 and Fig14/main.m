close all;
clear;
clc;

%% ------------------------------------------------------
% Provide the filename
% ------------------------------------------------------
% Radar raw data captured by TI AWR1243BOOST + TI mmwave DevPack + TI TSW1400
% TI AWR1243BOOST:      http://www.ti.com/tool/AWR1243BOOST
% TI mmwave DevPack:    http://www.ti.com/tool/MMWAVE-DEVPACK  
% TI TSW1400:           http://www.ti.com/tool/TSW1400EVM
% Configuration:        https://training.ti.com/sites/default/files/docs/mmWave_sensor_raw_data_capture_using_TSW1400_board_v3.pdf
fname = '20180608-195044.bin'; % With interference
fname2 = '20180613-190347.bin'; % Without interference
% fname2 = '20180613-190350.bin'; % Without interference
% fname2 = '../../.././Experiment/Bin/20180613-190353.bin';

%% ------------------------------------------------------
% Waveform Parameters
% ------------------------------------------------------
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

%% Open file of data with interference
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

%% Open file of data without interference
fid     = fopen(fname2,'r');     % Open the file in read mode
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
adc_dat2 = RX1;

%% Plot 10th chirp
OneChirp                = adc_dat(:, 10);  % Select the 50th chirp
OneRangeFFT             = fft(OneChirp, NFast);
priCH                   = OneRangeFFT(1:NFast/2);
refCH                   = conj(flipud(OneRangeFFT(NFast/2+1:NFast)));
RefPower                = sum(abs(refCH.^2))/(NFast/2);
%% SNR calculation
TargetIdx       = 74;
RefCells        = 20;
GuardCells      = 6;
SIRValue    = SIRCacl(priCH, TargetIdx, RefCells, GuardCells);
disp('*********************************************');
disp(sprintf('SIR in 10th chirp: %s', num2str(SIRValue)));
figure
set(gcf, 'Position', [0 0 1000 800])
subplot(421)
plot(real(OneChirp));
xlim([1 NFast]);
ylim([-0.03 0.03]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(a)');
xlabel(['Time (' 956 's)']);
ylabel('Amplitude(linear)');
xticklabels({'5', '10', '15', '20', '25'});
subplot(423)
plot(db(abs(priCH)));
xlim([1 240]);
ylim([-30 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(c)');
xlabel('Range(m)');
xticklabels({'10', '20', '30', '40'});
ylabel('Amplitude(dB)');
subplot(425)
plot(db(abs(refCH)));
xlim([1 240]);
ylim([-30 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(e)');
xlabel('Range(m)');
xticklabels({'10', '20', '30', '40'});
ylabel('Amplitude(dB)');


%% Setp size fraction
gamma = 20;
flen = 8;

%% Process 50th chirp
OneChirp                = adc_dat(:, 50);  % Select the 50th chirp
delay                   = -8;
OneChirp                = circshift(OneChirp,-delay);
OneRangeFFT             = fft(OneChirp, NFast);
priCH                   = OneRangeFFT(1:NFast/2);
refCH                   = conj(flipud(OneRangeFFT(NFast/2+1:NFast)));

priCH                   = flip(priCH);
% refCH                   = circshift(refCH,1);
refCH                   = flip(refCH);

RefPower                = sum(abs(refCH.^2))/(NFast/2);
% Adaptive filter
M                       = flen;                 % Length of adaptive filter
mu                      = 2/RefPower/gamma;    % Step size
[e1, wo]                = af_with_ref_cx(priCH, refCH, M, mu);

priCH                   = flip(priCH);
refCH                   = flip(refCH);
e1                      = flip(e1);

%% SNR calculation
TargetIdx       = 74;
RefCells        = 20;
GuardCells      = 6;
SIRbeforeANC    = SIRCacl(priCH, TargetIdx, RefCells, GuardCells);
disp('*********************************************');
disp(sprintf('SIR in 50th chirp before ANC: %s', num2str(SIRbeforeANC)));
SIRafterANC    = SIRCacl(e1, TargetIdx+1, RefCells, GuardCells);
disp(sprintf('SIR in 50th chirp after ANC: %s', num2str(SIRafterANC)));
disp(sprintf('SIR increased by (dB) in 50th chirp: %s \0', num2str(SIRafterANC - SIRbeforeANC)));

%% Plot
set(gcf, 'Position', [0 0 1000 800])
subplot(422)
plot(real(OneChirp));
xlim([1 NFast]);
ylim([-0.5 0.5]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(b)');
xlabel(['Time (' 956 's)']);
ylabel('Amplitude(linear)');
xticklabels({'5', '10', '15', '20', '25'});
subplot(424)
plot(db(abs(priCH)));
xlim([1 240]);
ylim([-30 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(d)');
xlabel('Range(m)');
xticklabels({'10', '20', '30', '40'});
ylabel('Amplitude(dB)');
subplot(426)
plot(db(abs(refCH)));
xlim([1 240]);
ylim([-30 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(f)');
xlabel('Range(m)');
xticklabels({'10', '20', '30', '40'});
ylabel('Amplitude(dB)');
subplot(428);
plot(db(abs(e1)));
xlim([1 240]);
ylim([-30 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(g)');
xlabel('Range(m)');
% xticklabels({'4', '8', '12', '16', '20', '24', '28', '32', '36', '40', '44', '48'});
xticklabels({'10', '20', '30', '40'});
ylabel('Amplitude(dB)');

% figure
% subplot(411)
% plot(angle(priCH));
% subplot(412)
% plot(angle((refCH)));
% subplot(413)
% plot(angle(priCH)-angle((refCH)));
% subplot(414)
% plot(abs(priCH-(refCH)));

% %% Time domain plots: just for test
% figure
% subplot(311);
% plot(real(OneChirp));        % Plot the received signal in the time domain
% xlim([1 NFast]);
% ylim([-0.5 0.5]);
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'bold';
% title('(a)');
% xlabel(['Time (' 956 's)']);
% ylabel('Amplitude(linear)');
% xticklabels({'5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});
% subplot(312);
% plot(real(ifft(priCH, NFast)));        % Plot the received signal in the time domain
% xlim([1 NFast]);
% ylim([-0.5 0.5]);
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'bold';
% title('(b)');
% xlabel(['Time (' 956 's)']);
% ylabel('Amplitude(linear)');
% xticklabels({'5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});
% subplot(313);
% plot(real(ifft(e1, NFast)));        % Plot the received signal in the time domain
% xlim([1 NFast]);
% ylim([-0.5 0.5]);
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'bold';
% title('(c)');
% xlabel(['Time (' 956 's)']);
% ylabel('Amplitude(linear)');
% xticklabels({'5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});


%% Plot Rang-Doppler before Adaptive Noise Canceler Without MTI
RangeFFT            = fft(adc_dat,NFast,1);                    % Fast time FFT
RDBeforeANC         = fftshift((fft(RangeFFT, NPRI, 2)));      % Slow time FFT
RDBeforeANC         = RDBeforeANC(NFast/2+1:NFast, :);            % Only with the positive RangeFFT part

%% Plot Rang-Doppler after Adaptive Noise Canceler without MTI
threshold = 0.5; % If the reference input power is greater than the threshold, there must be interferences.
for i = 1 : NPRI
    priCH                     = RangeFFT(1:NFast/2, i);
    refCH                     = conj(flipud( RangeFFT(NFast/2+1 : NFast , i)));  
    
    priCH               = flip(priCH);
    % refCH               = circshift(refCH,1);
    refCH               = flip(refCH);
    
    RefPower                  = sum(abs(refCH).^2)/(NFast/2);
    if (RefPower > threshold)
        % With adaptive filter
        M   = flen;      % Length of adaptive filter
        mu  = 2/RefPower/50;   % Weight update constant
        [e, wo] = af_with_ref_cx(priCH, refCH, M, mu);
        e                  =flip(e);
        eout(:, i) = e;
    else
        % Without Apaptive filter
        priCH               =flip(priCH);
        eout(:, i) = priCH;
    end
end
RDAfterANC          = fftshift(fft(eout, NPRI, 2), 2);
    
%% SNR Calculation
TargetRangeIdx      = 74;   % Get from the figure
TargetDopplerIdx    = 43;   % Get from the figure
RefCells            = 20;
GuardCells          = 6;
% Calculate SNR before ANC
NoisePwr            = 0;
for i = (TargetRangeIdx-RefCells/2-GuardCells/2) : (TargetRangeIdx+RefCells/2+GuardCells/2)
    for j = (TargetDopplerIdx-RefCells/2-GuardCells/2) : (TargetDopplerIdx+RefCells/2+GuardCells/2)
        if( (i<(TargetRangeIdx-GuardCells/2)||i>(TargetRangeIdx+GuardCells/2)) && ...
                (j<(TargetDopplerIdx-GuardCells/2)||j>(TargetDopplerIdx+GuardCells/2)) )
            NoisePwr = NoisePwr + abs(RDBeforeANC(i, j))^2;
        end
    end
end
NoisePwr = NoisePwr/((RefCells+1)^2-(GuardCells+1)^2);
TargetPwr           = abs(RDBeforeANC(TargetRangeIdx, TargetDopplerIdx))^2;
SIRBeforeANC        = 10*log10(TargetPwr/NoisePwr);
disp('*********************************************');
disp(sprintf('SIR before ANC in range-Doppler: %s \0', num2str(SIRBeforeANC)));
% Calculate SNR after ANC
NoisePwr            = 0;
for i = (TargetRangeIdx-RefCells/2-GuardCells/2) : (TargetRangeIdx+RefCells/2+GuardCells/2)
    for j = (TargetDopplerIdx-RefCells/2-GuardCells/2) : (TargetDopplerIdx+RefCells/2+GuardCells/2)
        if( (i<(TargetRangeIdx-GuardCells/2)||i>(TargetRangeIdx+GuardCells/2)) && ...
                (j<(TargetDopplerIdx-GuardCells/2)||j>(TargetDopplerIdx+GuardCells/2)) )
            NoisePwr = NoisePwr + abs(RDAfterANC(i, j))^2;
        end
    end
end
NoisePwr = NoisePwr/((RefCells+1)^2-(GuardCells+1)^2);
TargetPwr           = abs(RDAfterANC(TargetRangeIdx, TargetDopplerIdx))^2;
SIRAfterANC        = 10*log10(TargetPwr/NoisePwr);
disp(sprintf('SIR after ANC in range-Doppler: %s \0', num2str(SIRAfterANC)));
% Calculate increased SNR
disp(sprintf('SIR increased by (dB) without MTI: %s \0', num2str(SIRAfterANC - SIRBeforeANC)));

%%
zlimup = 40;
zlimdown = -10;

%% Scale to range and velocity
MaxRange                = 46; % 40 meter
MaxRangeFFTIdx          = round(2*MaxRange/c*ChirpRate*NFast/ADCFs);
MaxVelocity             = 50; % MPH
MaxVelocity             = MaxVelocity*0.44704; % m/s
MaxVelocityIdx          = round(2*MaxVelocity*Fc/c*NPRI*PRI);
[DopplerIdx, RangeIdx]  = meshgrid(1:MaxVelocityIdx*2+1, 1:MaxRangeFFTIdx);
RangeIdx                = RangeIdx/(NFast/2)*(ADCFs/2)/ChirpRate*c/2;
RangeScaled             = [1 MaxRangeFFTIdx];
RangeScaled             = RangeScaled/(NFast/2)*(ADCFs/2)/ChirpRate*c/2;
DopplerIdx              = (DopplerIdx - MaxVelocityIdx -1)/NPRI/PRI*c/Fc/2/0.44704;
VelocityScaled          = [-MaxVelocityIdx MaxVelocityIdx-1];
VelocityScaled          = VelocityScaled/NPRI/PRI*c/Fc/2/0.44704; % mph

% Plot
RangeFFT2               = fft(adc_dat2,NFast,1);                    % Fast time FFT
RDBeforeANC2            = fftshift((fft(RangeFFT2, NPRI, 2)));      % Slow time FFT
RDBeforeANC2            = RDBeforeANC2(NFast/2+1:NFast, :);            % Only with the positive RangeFFT part
RangeDopplerScaled2     = RDBeforeANC2(1 : MaxRangeFFTIdx, NPRI/2+1-MaxVelocityIdx: NPRI/2+1+MaxVelocityIdx);
TargetRangeIdx      = 72;   % Get from the figure
TargetDopplerIdx    = 40;   % Get from the figure
RefCells            = 20;
GuardCells          = 6;
% Calculate SNR before ANC
NoisePwr            = 0;
for i = (TargetRangeIdx-RefCells/2-GuardCells/2) : (TargetRangeIdx+RefCells/2+GuardCells/2)
    for j = (TargetDopplerIdx-RefCells/2-GuardCells/2) : (TargetDopplerIdx+RefCells/2+GuardCells/2)
        if( (i<(TargetRangeIdx-GuardCells/2)||i>(TargetRangeIdx+GuardCells/2)) && ...
                (j<(TargetDopplerIdx-GuardCells/2)||j>(TargetDopplerIdx+GuardCells/2)) )
            NoisePwr = NoisePwr + abs(RDBeforeANC2(i, j))^2;
        end
    end
end
NoisePwr = NoisePwr/((RefCells+1)^2-(GuardCells+1)^2);
TargetPwr           = abs(RDBeforeANC2(TargetRangeIdx, TargetDopplerIdx))^2;
SIRBeforeANC2       = 10*log10(TargetPwr/NoisePwr);
disp('*********************************************');
disp(sprintf('SIR in range-Doppler without interference: %s \0', num2str(SIRBeforeANC2)));
disp('*********************************************');
dBRangeDopplerScaled2 = db(abs(RangeDopplerScaled2));
dBRangeDopplerScaled2(dBRangeDopplerScaled2<zlimdown) = zlimdown;
figure
% subplot(221);
surf(DopplerIdx, RangeIdx, dBRangeDopplerScaled2, 'LineStyle', 'none');
colorbar;
caxis([zlimdown zlimup]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
xlabel('Speed(mph)');
xlim(VelocityScaled);
ylabel('Range(meter)');
ylim(RangeScaled);
zlim([zlimdown zlimup]);
view(2);
view([46 68]);

% Plot
figure
% subplot(222)
RangeDopplerScaled      = RDBeforeANC(1 : MaxRangeFFTIdx, NPRI/2+1-MaxVelocityIdx: NPRI/2+1+MaxVelocityIdx);
dBRangeDopplerScaled = db(abs(RangeDopplerScaled));
dBRangeDopplerScaled(dBRangeDopplerScaled<zlimdown) = zlimdown;
surf(DopplerIdx, RangeIdx, dBRangeDopplerScaled, 'LineStyle', 'none');
colorbar;
caxis([zlimdown zlimup]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
xlabel('Speed(mph)');
xlim(VelocityScaled);
ylabel('Range(meter)');
ylim(RangeScaled);
zlim([zlimdown zlimup]);
view(2);
view([60 68]);

figure
% subplot(212);
RangeDopplerScaled      = RDAfterANC(1 : MaxRangeFFTIdx, NPRI/2+1-MaxVelocityIdx: NPRI/2+1+MaxVelocityIdx);
dBRangeDopplerScaled = db(abs(RangeDopplerScaled));
dBRangeDopplerScaled(dBRangeDopplerScaled<zlimdown) = zlimdown;
surf(DopplerIdx, RangeIdx, dBRangeDopplerScaled, 'LineStyle', 'none');
colorbar;
caxis([zlimdown zlimup]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
xlabel('Speed(mph)');
xlim(VelocityScaled);
ylabel('Range(meter)');
ylim(RangeScaled);
zlim([zlimdown zlimup]);
view(2);
view([53 71]);