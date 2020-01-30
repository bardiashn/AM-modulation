% In the name of wisdom
%==========================================================================
% Student Name: amirhossein shoaraye nejati
% ID: 9521010039
% Project 2
% File name: Am_Modulator.m(by a Sawtooth function input)
% Date: 1.10.2019
% Kashan University
% Course: Fundamentals of Communications Systems, Fall 2018
%==========================================================================
%%%%%%%parameter%%%%%%
clc;% Clear Command Window
clear;    % Remove items from workspace, freeing up system memory
close all;    % close all figures
SNR=0.01; %signal to noise ratio
fm=10e2;%massage frequancey 
fc=2e4;% carrier frequancey
Am=20000;%massage amplitude volt
Ac=10;% carrier amplitude volt
FS=1e6;     % Fs=1 MHz %Sampling Frequency (Hz) % Should be at least 2 x (highest freq. of all signals)
TS=1/FS;      % Sampling Time
n=2^20;       % number of fft points 
g=1/Ac^2;      %gain of filter
B=1;            % AM Mod. index
f=-FS/2:FS/n:(FS/2-FS/n);% frequency variable
%%%%%%% product modulator%%%%%%
t=-0.1:TS:0.1;  % time variable
massage = sawtooth(2*pi*100*t,0.2); % Sawtooth signal
m_f = TS*fftshift(fft(massage,n));
%carrier%
carrier=Ac*cos(2*pi*fc*t);
c_f = TS*fftshift(fft(carrier,n));
%transmitter signal
AM_signal=(1+B*massage).*carrier;
a_f= TS*fftshift(fft(AM_signal,n));
%%%%%%% channel%%%%%%%
L=1; % Channel power loss 
Sigma_n2=20 ; % noise variance (set Sigma_n2=0 for No noise)
n_t=sqrt(Sigma_n2) * randn(1,length(AM_signal)); % n(t) = AWGN noise
tX_signal= 1/sqrt(L) * AM_signal  + n_t ;  % yR(t) = Received signal at receiver

%%%%%%product demodulater%%%%%%
%recevied signal%
RX_signal =  tX_signal.*carrier;
r_f = TS*fftshift(fft(RX_signal,n));
% low pass filter%
wn = .02;     % PSD lowpass filter cut - off frequency
[b,a] = butter(2,wn);    % Design lowpass filter
demod_signal =2*((filter(b,a,RX_signal).*g)-0.5);     % Apply lowpass filter
d_f = TS*fftshift(fft(demod_signal,n));                                             
%%%%PLOTS%%%%%
%plot massage signal%
figure
 subplot(2,1,1)
 plot(t,massage,'r')
 grid on
 title('massage signal')
 xlabel('time(s)')
 ylabel('amplitude ')
 ylim([-1.1 1.1])
 legend ('time domain ---->');
 subplot(2,1,2)
 plot(f,abs(m_f),'r')
 title('massage signal in frequency domain')
 xlabel('frequency(Hz)')
 ylabel('magnitude')
 xlim([-1000 1000]);
 ylim([0 0.1]);
 legend ('frequency domain---->');
 %plot carrier signal%
 figure
 subplot(2,1,1)
 plot(t,carrier,'g')
 grid on
 title('carrier signal')
 xlabel('time(s)')
 ylabel('amplitude ')
 xlim([-1.5e-3 1.5e-3])
 ylim([7 10])
 legend('time domain ---->');
 subplot(2,1,2)
 plot(f,abs(c_f),'g')
 title('carrier signal in frequency domain')
 xlabel('frequency(Hz)')
 ylabel('magnitude')
 ylim([0 1.2])
 xlim([-0.4e5 0.4e5])
 legend('freuency domain---->');
 %plot transmitted signal%
 figure
 subplot(2,1,1)
  plot(t,AM_signal,'b')
 grid on
 title('transmitted signal')
 xlabel('time(s)')
 ylabel('amplitude ')
 legend('time domain ---->');
 subplot(2,1,2)
 plot(f,abs(a_f),'b')
 title('transmitted signal in frequency domain')
 xlabel('frequency(hz)')
 ylabel('magnitude')
 ylim([0 1.2])
 xlim([-3e4 3e4])
 legend('frequency domain --->')
 %plot recevied signal%
 figure
 subplot(2,1,1)
  plot(t,RX_signal,'k')
 grid on
 title('recevied signal')
 xlabel('time(s)')
 ylabel('amplitude ')
 ylim([0 250])
 xlim([-0.12 0.12])
 legend('time domain---->');
 subplot(2,1,2)
 plot(f,abs(r_f),'k')
 title('recevied signal in frequency domain')
 xlabel('frequency(hz)')
 ylabel('magnitude')
 ylim([0 15])
 xlim([-500 500])
 legend('frequency domain ---->');
 %plot demodulated signl%
  figure
 subplot(2,1,1)
  plot(t,demod_signal,'m')
 grid on
 title('Demodulated  signal')
 xlabel('time(s)')
 ylabel('amplitude ')
 ylim([-1.1 1.1])
 xlim([-0.05 0.05])
 legend('time domain---->');
 subplot(2,1,2)
 plot(f,abs(d_f),'m')
 title('demodulated signal in frequency domain')
 xlabel('frequency(Hz)')
 ylabel('magnitude')
 ylim ([0 0.1])
 xlim([-1e4 1e4])
 legend('freuency domain---->');
 
%END
%==========================================================================

 