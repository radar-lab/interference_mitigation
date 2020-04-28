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
% fname = '20180613-190350.bin'; % Without interference
fname2 = '20180613-190350.bin'; % Without interference
% fname2 = '../../.././Experiment/Bin/20180613-190325.bin';

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

%% Open file and pre-porcess the data
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
ylim([-40 20]);
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
ylim([-40 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(e)');
xlabel('Range(m)');
xticklabels({'10', '20', '30', '40'});
ylabel('Amplitude(dB)');


%% Setp size fraction
gamma = 30;
flen = 8;

%% Process 50th chirp
OneChirp                = adc_dat(:, 50);  % Select the 50th chirp
OneRangeFFT             = fft(OneChirp, NFast);
priCH                   = OneRangeFFT(1:NFast/2);
refCH                   = conj(flipud(OneRangeFFT(NFast/2+1:NFast)));
RefPower                = sum(abs(refCH.^2))/(NFast/2);
% Adaptive filter
M                       = flen;                 % Length of adaptive filter
mu                      = 2/RefPower/gamma;    % Step size
[e1, wo]                = af_with_ref_cx(priCH, refCH, M, mu);

%% SNR calculation
TargetIdx       = 74;
RefCells        = 20;
GuardCells      = 6;
increasedSNR    = IncreasedSNR(priCH, e1, TargetIdx, RefCells, GuardCells);
disp(sprintf('SNR increased by (dB) for mu = 30: %s \0', num2str(increasedSNR)));

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
ylim([-40 20]);
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
ylim([-40 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(f)');
xlabel('Range(m)');
xticklabels({'10', '20', '30', '40'});
ylabel('Amplitude(dB)');
subplot(414);
plot(db(abs(e1)));
xlim([1 240]);
ylim([-40 20]);
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
title('(g)');
xlabel('Range(m)');
xticklabels({'4', '8', '12', '16', '20', '24', '28', '32', '36', '40', '44', '48'});
ylabel('Amplitude(dB)');


%% Plot Rang-Doppler before Adaptive Noise Canceler Without MTI
RangeFFT            = fft(adc_dat,NFast,1);                    % Fast time FFT
RDBeforeANC         = fftshift((fft(RangeFFT, NPRI, 2)));      % Slow time FFT
RDBeforeANC         = RDBeforeANC(NFast/2+1:NFast, :);            % Only with the positive RangeFFT part

%% Plot Rang-Doppler before Adaptive Noise Canceler With MTI
h = [1 -1];
RangeFFTMTI         = fft(adc_dat,NFast,1);                    % Fast time FFT
RangeFFTMTI         = filter(h,1,RangeFFTMTI,[],2);
RDBeforeANCMTI      = fftshift((fft(RangeFFTMTI, NPRI, 2)));      % Slow time FFT
RDBeforeANCMTI      = RDBeforeANCMTI(NFast/2+1:NFast, :);            % Only with the positive RangeFFT part

%% Plot Rang-Doppler after Adaptive Noise Canceler without MTI
threshold = 0.5; % If the reference input power is greater than the threshold, there must be interferences.
for i = 1 : NPRI
    priCH                     = RangeFFT(1:NFast/2, i);
    refCH                     = conj(flipud( RangeFFT(NFast/2+1 : NFast , i)));   
    RefPower                  = sum(abs(refCH).^2)/(NFast/2);
    if (RefPower > threshold)
        % With adaptive filter
        M   = flen;      % Length of adaptive filter
        mu  = 2/RefPower/gamma;   % Weight update constant
        [e, wo] = af_with_ref_cx(priCH, refCH, M, mu);
        eout(:, i) = e;
    else
        % Without Apaptive filter
        eout(:, i) = priCH;
    end
end
RDAfterANC          = fftshift(fft(eout, NPRI, 2), 2);

%% Plot Rang-Doppler after Adaptive Noise Canceler with MTI
threshold = 0.5; % If the reference input power is greater than the threshold, there must be interferences.
for i = 1 : NPRI
    priCH                     = RangeFFTMTI(1:NFast/2, i);
    refCH                     = conj(flipud( RangeFFTMTI(NFast/2+1 : NFast , i)));   
    RefPower                  = sum(abs(refCH).^2)/(NFast/2);
    if (RefPower > threshold)
        % With adaptive filter
        M   = flen;      % Length of adaptive filter
        mu  = 2/RefPower/gamma;   % Weight update constant
        [e, wo] = af_with_ref_cx(priCH, refCH, M, mu);
        eout(:, i) = e;
    else
        % Without Apaptive filter
        eout(:, i) = priCH;
    end
end
RDAfterANCMTI          = fftshift(fft(eout, NPRI, 2), 2);
    
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
% Calculate SNR after ANC without MTI
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
% Calculate increased SNR
disp(sprintf('SNR increased by (dB) without MTI:: %s \0', num2str(SIRAfterANC - SIRBeforeANC)));
% Calculate SNR after ANC with MTI
NoisePwr            = 0;
for i = (TargetRangeIdx-RefCells/2-GuardCells/2) : (TargetRangeIdx+RefCells/2+GuardCells/2)
    for j = (TargetDopplerIdx-RefCells/2-GuardCells/2) : (TargetDopplerIdx+RefCells/2+GuardCells/2)
        if( (i<(TargetRangeIdx-GuardCells/2)||i>(TargetRangeIdx+GuardCells/2)) && ...
                (j<(TargetDopplerIdx-GuardCells/2)||j>(TargetDopplerIdx+GuardCells/2)) )
            NoisePwr = NoisePwr + abs(RDAfterANCMTI(i, j))^2;
        end
    end
end
NoisePwr = NoisePwr/((RefCells+1)^2-(GuardCells+1)^2);
TargetPwr           = abs(RDAfterANCMTI(TargetRangeIdx, TargetDopplerIdx))^2;
SIRAfterANCMTI        = 10*log10(TargetPwr/NoisePwr);
% Calculate increased SNR
disp(sprintf('SNR increased by (dB) with MTI:: %s \0', num2str(SIRAfterANCMTI - SIRBeforeANC)));


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
figure
% subplot(221)
RangeFFT2            = fft(adc_dat2,NFast,1);                    % Fast time FFT
RDBeforeANC2         = fftshift((fft(RangeFFT2, NPRI, 2)));      % Slow time FFT
RDBeforeANC2         = RDBeforeANC2(NFast/2+1:NFast, :);            % Only with the positive RangeFFT part
RangeDopplerScaled2      = RDBeforeANC2(1 : MaxRangeFFTIdx, NPRI/2+1-MaxVelocityIdx: NPRI/2+1+MaxVelocityIdx);
dBRangeDopplerScaled2 = db(abs(RangeDopplerScaled2));
dBRangeDopplerScaled2(dBRangeDopplerScaled2<zlimdown) = zlimdown;
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

% Plot
figure
% subplot(221)
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
figure
% subplot(222);
RangeDopplerScaled      = RDBeforeANCMTI(1 : MaxRangeFFTIdx, NPRI/2+1-MaxVelocityIdx: NPRI/2+1+MaxVelocityIdx);
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
figure
% subplot(223);
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
figure
% subplot(224);
RangeDopplerScaled      = RDAfterANCMTI(1 : MaxRangeFFTIdx, NPRI/2+1-MaxVelocityIdx: NPRI/2+1+MaxVelocityIdx);
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