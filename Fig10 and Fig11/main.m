close all;
clear;
clc;

%% Common parameters
SystemFs        = 1e9;                      % Simulation resolution in the time domain
c               = 3e8;                      % Speed of light

%% Waveform parameters
Fc              = 76e9;                     % Carrier Frequency
BW              = 300e6;                    % Bandwidth of chirp signal
PRI             = 51.2e-6;                  % Chirp Duration
ChirpRate       = BW/PRI;                   % Chirp Rate
ADCFs           = 40e6;                     % ADC sampling frequency
N_Fast          = round(PRI*ADCFs);         % Number of the ADC samples during one chirp
FFTSize_Fast    = 2^(round(log2(N_Fast)));  % Number of the fast time FFT size
PRISize         = PRI*SystemFs;             % Number of the data size of one chirp

%% Time scale
t  = (0 : 1/SystemFs : PRI-1/SystemFs);     % Generate the time instant

%% TX and ref signal
TX_RF   = exp(1i*pi*ChirpRate*t.^2).*exp(1i*2*pi*Fc*t); % Generate one transmitting chirp signal
TX_Ref  = conj(TX_RF);                                  % Genertae the LO signal for the mixer

%% Echo #1 generation
T1_R    = 40;       % Target at distance of 25 meters
Echo1   = 2*exp(1i*pi*ChirpRate*(t-2*T1_R/c).^2).*exp(1i*2*pi*Fc*(t-2*T1_R/c)) / ((2*T1_R)^2); % Target echo with amplitude attunation effect due to the distance

%% Echo #2 generation
T2_R    = 100;       % Target at distance of 100 meters
Echo2   = 5*exp(1i*pi*ChirpRate*(t-2*T2_R/c).^2).*exp(1i*2*pi*Fc*(t-2*T2_R/c)) / ((2*T2_R)^2); % Target echo with amplitude attunation effect due to the distance

%% Interference #1 generation
IF1_R       = 10;   % Interference #1 at distance of 10 meters
FactorN     = 4;    % Interference #1's chirp rate is 4 times of the transmitting chirp
IF1_CR      = FactorN * ChirpRate;
for i = 1 : FactorN
    t_sec = (i-1)*PRI/FactorN  : 1/SystemFs : i*PRI/FactorN-1/SystemFs;
    IF1_Sig(1+(i-1)*PRISize/FactorN : i*PRISize/FactorN) = exp(1i*2*pi*(-BW*(i-1))*t_sec).*exp(1i*pi*IF1_CR*t_sec.^2).*exp(1i*2*pi*Fc*t_sec);
end
IF1_Sig     = IF1_Sig / (IF1_R^2); % With amplitude attunation effect due to the distance

%% Interference #2 generation
IF2_R       = 30;   % Interference #2 at distance of 20 meters
FactorN     = 5;    % Interference #2's chirp rate is 5 times of the transmitting chirp
IF2_CR      = FactorN * ChirpRate;
for i = 1 : FactorN
    t_sec = (i-1)*PRI/FactorN : 1/SystemFs : i*PRI/FactorN-1/SystemFs;
    IF2_Sig(1+(i-1)*PRISize/FactorN : i*PRISize/FactorN) = exp(1i*2*pi*(-BW*(i-1))*t_sec).*exp(1i*pi*IF2_CR*t_sec.^2).*exp(1i*2*pi*Fc*t_sec);
end
IF2_Sig     = IF2_Sig / (IF2_R^2); % With amplitude attunation effect due to the distance

%% Interfernece #3 generation
IF3_R       = 50;   % Interference #3 at distance of 30 meters
IF3_Sig     = exp(1i*2*pi*100e6*t).*exp(1i*2*pi*Fc*t) / (IF3_R^2); % CW interference with amplitude attunation effect due to the distance

%% Received RF signal generation
RX_Sig      = Echo1 + Echo2 + IF1_Sig + IF2_Sig + IF3_Sig; % Received signal consistes of one echo and three interferences

%% LNA
Gain        = 40; % in dB
RX_Sig      = RX_Sig * 10^(Gain/20);

%% Mixer
Mixer_Output = conj(RX_Sig.*TX_Ref);

%% Analog LPF
LPF_Output  = filter(FIR, Mixer_Output); % Pass band is 10MHz, Stop band is 20MHz

%% ADC
RX_Base     = downsample(LPF_Output, round(SystemFs/ADCFs)); % Use downsampling to include the ADC effect

%% Adaptive Noise Canceler (ANC)
RangeDoppler        = fftshift((fft(RX_Base,FFTSize_Fast)));            % FFT of the received baseband signal RX_Base
priCH               = RangeDoppler(FFTSize_Fast/2+1 : FFTSize_Fast);    % Generate the primary channel input for the ANC
refCH               = conj(fliplr(RangeDoppler(1 : FFTSize_Fast/2)));   % Generate the reference channel input for the ANC
priCH               = priCH.';                                          % Make it M-by-1 vector
refCH               = refCH.';                                          % Make it M-by-1 vector
RefPower            = sum(abs(refCH).^2)/(FFTSize_Fast/2);              % Estimate the input power of the reference channel
M                   = 8;                                                % Length of adaptive filter
mu                  = 2/RefPower/50;                                    % Step size mu
[e1, wo1]           = af_with_ref_cx(priCH, refCH, M, mu);              % Apply complex ANC, e is the final output
mu                  = 2/RefPower/100;                                    % Step size mu
[e2, wo2]           = af_with_ref_cx(priCH, refCH, M, mu);              % Apply complex ANC, e is the final output
mu                  = 2/RefPower/150;                                    % Step size mu
[e3, wo3]           = af_with_ref_cx(priCH, refCH, M, mu);              % Apply complex ANC, e is the final output

%% Plot the STFT of the received signal before and after the mixer
STFT_WindonwLen = 512;
figure
set(gcf, 'Position', [0 0 1000 800])
subplot(311)
% Plot the STFT of the received signal before the mixer
spectrogram(RX_Sig,STFT_WindonwLen,round(STFT_WindonwLen*0.8),STFT_WindonwLen, SystemFs, 'centered','yaxis');
ylim([-100 400]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
title('(a)');
STFT_WindonwLen = 1024;
subplot(312)
% Plot the STFT of the received signal after the mixer
spectrogram(LPF_Output,STFT_WindonwLen,round(STFT_WindonwLen*0.8),STFT_WindonwLen, SystemFs, 'centered','yaxis');
ylim([-20 20]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
title('(b)');
subplot(313);
plot(real(RX_Base));        % Plot the received signal in the time domain
xlim([1 N_Fast]);
ylim([-1.5 1.5]);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
title('(c)');
xlabel(['Time (' 956 's)']);
ylabel('Amplitude(linear)');
xticklabels({'5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});

% %% Time domian plots: just for test
% figure
% subplot(311);
% plot(real(RX_Base));        % Plot the received signal in the time domain
% xlim([1 N_Fast]);
% ylim([-1.5 1.5]);
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'bold';
% title('(a)');
% xlabel(['Time (' 956 's)']);
% ylabel('Amplitude(linear)');
% xticklabels({'5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});
% subplot(312);
% plot(real(ifft(priCH, FFTSize_Fast)));        % Plot the received signal in the time domain
% xlim([1 N_Fast]);
% ylim([-1.5 1.5]);
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'bold';
% title('(b)');
% xlabel(['Time (' 956 's)']);
% ylabel('Amplitude(linear)');
% xticklabels({'5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});
% subplot(313);
% plot(real(ifft(e1, FFTSize_Fast)));        % Plot the received signal in the time domain
% xlim([1 N_Fast]);
% ylim([-1.5 1.5]);
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'bold';
% title('(c)');
% xlabel(['Time (' 956 's)']);
% ylabel('Amplitude(linear)');
% xticklabels({'5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});

%% Plot the simulation result of the adaptive noise canceler
figure
set(gcf, 'Position', [0 0 1000 800])
subplot(511);
plot(db(abs(priCH)));         % Plot the positive FFT part of the received signal
xlim([1 600]);
xticklabels({'25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '275', '300'});
ylim([0 40]);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('(a)');
xlabel('Range(m)');
ylabel('Amplitude(dB)');
subplot(512);
plot(db(abs(refCH)));         % Plot the flipped conjugated negative FFT part of the received signal
xlim([1 600]);
xticklabels({'25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '275', '300'});
ylim([0 40]);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('(b)');
xlabel('Range(m)');
ylabel('Amplitude(dB)');
% figure
subplot(513);
plot(db(abs(e1)));             % Plot the output of the ANC
xlim([1 600]);
xticklabels({'25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '275', '300'});
ylim([0 40]);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('(c)');
xlabel('Range(m)');
ylabel('Amplitude(dB)');
subplot(514);
plot(db(abs(e2)));             % Plot the output of the ANC
xlim([1 600]);
xticklabels({'25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '275', '300'});
ylim([0 40]);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('(d)');
xlabel('Range(m)');
ylabel('Amplitude(dB)');
subplot(515);
plot(db(abs(e3)));             % Plot the output of the ANC
xlim([1 600]);
xticklabels({'25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '275', '300'});
ylim([0 40]);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
title('(e)');
xlabel('Range(m)');
ylabel('Amplitude(dB)');

%% SNR calculation
TargetIdx = T1_R/(c/2/BW)+1; % Target index
RefCells    = 20;
GuardCells   = 6;
disp('Target #1:');
IncreasedSNR1 = IncreasedSNR(priCH, e1, TargetIdx, RefCells, GuardCells);
disp(sprintf('SNR increased by (dB) for mu = 50: %s \0', num2str(IncreasedSNR1)));
IncreasedSNR2 = IncreasedSNR(priCH, e2, TargetIdx, RefCells, GuardCells);
disp(sprintf('SNR increased by (dB) for mu = 100: %s \0', num2str(IncreasedSNR2)));
IncreasedSNR3 = IncreasedSNR(priCH, e3, TargetIdx, RefCells, GuardCells);
disp(sprintf('SNR increased by (dB) for mu = 150: %s \0', num2str(IncreasedSNR3)));

TargetIdx = T2_R/(c/2/BW)+1; % Target index
RefCells    = 20;
GuardCells   = 6;
disp('Target #2:');
IncreasedSNR1 = IncreasedSNR(priCH, e1, TargetIdx, RefCells, GuardCells);
disp(sprintf('SNR increased by (dB) for mu = 50: %s \0', num2str(IncreasedSNR1)));
IncreasedSNR2 = IncreasedSNR(priCH, e2, TargetIdx, RefCells, GuardCells);
disp(sprintf('SNR increased by (dB) for mu = 100: %s \0', num2str(IncreasedSNR2)));
IncreasedSNR3 = IncreasedSNR(priCH, e3, TargetIdx, RefCells, GuardCells);
disp(sprintf('SNR increased by (dB) for mu = 150: %s \0', num2str(IncreasedSNR3)));


% %% SNR vs mu plot
% SNRvsMu = zeros(1, 200);
% for i = 15 : 200
%     M                   = 8;                                                % Length of adaptive filter
%     mu                  = 2/RefPower/i;                                    % Step size mu
%     [e, wo]             = af_with_ref_cx(priCH, refCH, M, mu);              % Apply complex ANC, e is the final output
%     
%     RefCells    = 20;
%     GuardCells  = 6;
%     SNRvsMu(i)  = IncreasedSNR(priCH, e, TargetIdx, RefCells, GuardCells);
% end
% figure
% plot(SNRvsMu);
% xlabel('Fraction of the upper bound $\displaystyle\frac{2}{Filter Input Power}$','interpreter','latex', 'FontSize', 14, 'FontWeight','bold');
% ylabel('Increased SIR', 'FontSize', 18, 'FontWeight','bold');%, 'FontWeight','bold');
% set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'1', '$\frac{1}{20}$', '$\frac{1}{40}$', ...
%     '$\frac{1}{60}$', '$\frac{1}{80}$', '$\frac{1}{100}$', '$\frac{1}{120}$', '$\frac{1}{140}$', '$\frac{1}{160}$', '$\frac{1}{180}$', '$\frac{1}{200}$'}, 'FontWeight','bold');