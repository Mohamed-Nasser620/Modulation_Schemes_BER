%% Global Variables
SNR_Arr = [1 2 3 4 5 6 7 8 9 10]; %SNR Range
num_bits = 3600000; %number of bits in stream
 
%generating a stream of random bits
data_bits = randi([0,1], 1, num_bits);
 
%% BPSK Variables
Esymb_BPSK = 1; %symbol Energy = bit energy = 1
symbols_BPSK = zeros(1,num_bits); %array that will contain the symbols representing the bits
 
%generating normally distributed noise_BPSK with zero mean and variance = 1
noise_BPSK = randn(1, num_bits);
scaled_noise_BPSK = zeros(length(SNR_Arr), num_bits); %scaled noise_BPSK initialization
No_BPSK = zeros(1,length(SNR_Arr)); %No_BPSK initialization for each SNR
Rx_BPSK = zeros(length(SNR_Arr), num_bits); %Received data after noise_BPSK
%Received bits after comparing to descision regions using decision boundaries
Rx_BPSK_bits = zeros(length(SNR_Arr), num_bits);
%an array of vectors that identifies the error_BPSK occured in which bits for
%each received stream after adding noise to transmitted stream 
error_BPSK = zeros(length(SNR_Arr), num_bits); 
BER_BPSK = zeros(1,length(SNR_Arr)); %bit error rate for each stream of bits in BPSK
 
%% BPSK Mapping
%transfroming bits into symbols '1' to 1 volt and '0' to -1 volt
symbols_BPSK = (2*data_bits - 1) * Esymb_BPSK;
 
 
%% QPSK Variables
Esymb_QPSK = 2; %Symbol Energy = bit energy * 2
%array that will contain symbols representing bits using gray encoding
symbols_QPSK_GE = zeros(1,num_bits/2); 
%array that will contain symbols representing bits without using gray encoding
symbols_QPSK_NGE = zeros(1,num_bits/2);
 
%generating normally distributed complex noise_QPSK with zero mean and variance = 1
noise_QPSK = randn(1, num_bits/2)+1j*randn(1, num_bits/2);
scaled_noise_QPSK = zeros(length(SNR_Arr), num_bits/2); %scaled noise_QPSK initialization
No_QPSK = zeros(1,length(SNR_Arr)); %No_QPSK initialization for each SNR
 
%Gray Encoding Variables
Rx_QPSK_GE = zeros(length(SNR_Arr), num_bits/2); %Received data after noise_QPSK
%Received bits after comparing to descision regions using decision boundaries
Rx_QPSK_GE_bits = zeros(length(SNR_Arr), num_bits); 
%an array of vectors that identifies the error_QPSK_GE occured in which bits for
%each received stream after adding noise to transmitted stream 
error_QPSK_GE = zeros(length(SNR_Arr), num_bits); 
BER_QPSK_GE = zeros(1,length(SNR_Arr)); %bit error rate for each stream of bits
 
%Not Gray Encoding Variables
Rx_QPSK_NGE = zeros(length(SNR_Arr), num_bits/2); %Received data after noise_QPSK
%Received bits after comparing to descision regions using decision boundaries
Rx_QPSK_NGE_bits = zeros(length(SNR_Arr), num_bits); 
%an array of vectors that identifies the error_QPSK_NGE occured in which bits for
%each received stream after adding noise to transmitted stream 
error_QPSK_NGE = zeros(length(SNR_Arr), num_bits); 
BER_QPSK_NGE = zeros(1,length(SNR_Arr)); %bit error rate for each stream of bits
 
%% QPSK Mapping
%looping in stream of bits with step = 2 to represent each two bits with
%complex numbers to identify constellation points using gray encoding and
%non gray encoding
for n = 1 : 2 : num_bits
    if [data_bits(n) data_bits(n+1)] == [0 0]
        symbols_QPSK_GE(floor(n/2 +1))= -1-1j; % 00 represented by -1-j
        symbols_QPSK_NGE(floor(n/2 +1))= -1-1j;% 00 represented by -1-j
    elseif [data_bits(n) data_bits(n+1)] == [0 1]
        symbols_QPSK_GE(floor(n/2 +1))= -1+1j; % 01 represented by -1+j
        symbols_QPSK_NGE(floor(n/2 +1))= -1+1j;% 01 represented by -1+j
    elseif [data_bits(n) data_bits(n+1)] == [1 0]
        symbols_QPSK_GE(floor(n/2 +1))= 1-1j; % 10 represented by 1-j
        symbols_QPSK_NGE(floor(n/2 +1))= 1+1j;% 10 represented by 1+j
    elseif [data_bits(n) data_bits(n+1)] == [1 1]
        symbols_QPSK_GE(floor(n/2 +1))= 1+1j; % 11 represented by 1+j
        symbols_QPSK_NGE(floor(n/2 +1))= 1-1j;% 11 represented by 1-j
    end
end
 
 
%% 8PSK Variables
Esymb_8PSK = 1; %symbol Energy = bit energy * 3
symbols_8PSK = zeros(1,num_bits/3); %array that will contain the symbols representing the bits
 
%generating normally distributed complex noise_8PSK with zero mean and variance = 1
noise_8PSK = randn(1, num_bits/3)+1j*randn(1, num_bits/3);
scaled_noise_8PSK = zeros(length(SNR_Arr), num_bits/3); %scaled noise_8PSK initialization
No_8PSK = zeros(1,length(SNR_Arr)); %No_8PSK initialization for each SNR
Rx_8PSK = zeros(length(SNR_Arr), num_bits/3); %Received data after added noise_8PSK
%Received bits after comparing to descision regions using decision boundaries
Rx_8PSK_bits = zeros(length(SNR_Arr), num_bits); 
%an array of vectors that identifies the error_8PSK occured in which bits for
%each received stream after adding noise to transmitted stream 
error_8PSK = zeros(length(SNR_Arr), num_bits); 
BER_8PSK = zeros(1,length(SNR_Arr)); %bit error rate for each stream of bits
 
%% 8PSK Mapping
%looping in stream of bits with step = 3 to represent each 3 bits with
%complex numbers to identify constellation points using gray encoding 
for n = 1 : 3 : num_bits
    if [data_bits(n) data_bits(n+1) data_bits(n+2)] == [0 0 0]
        symbols_8PSK(floor(n/3 +1))= 1+0j; 
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2)] == [0 0 1]
        symbols_8PSK(floor(n/3 +1))= (sqrt(2)/2)+1j*(sqrt(2)/2);
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2)] == [0 1 0]
        symbols_8PSK(floor(n/3 +1))= -(sqrt(2)/2)+1j*(sqrt(2)/2);
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2)] == [0 1 1]
        symbols_8PSK(floor(n/3 +1))= 0+1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2)] == [1 0 0]
        symbols_8PSK(floor(n/3 +1))= (sqrt(2)/2)-1j*(sqrt(2)/2);
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2)] == [1 0 1]
        symbols_8PSK(floor(n/3 +1))= 0-1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2)] == [1 1 0]
        symbols_8PSK(floor(n/3 +1))= -1+0j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2)] == [1 1 1]
        symbols_8PSK(floor(n/3 +1))= -(sqrt(2)/2)-1j*(sqrt(2)/2);    
    end
end
 
 
%% 16QAM Variables
Esymb_16qam = 1; %symbol Energy = bit energy / 2.5
symbols_16qam = zeros(1,num_bits/4); %array that will contain the symbols representing the bits
 
%generating normally distributed complex noise_16qam with zero mean and variance = 1
noise_16qam = randn(1, num_bits/4)+1j*randn(1, num_bits/4);
scaled_noise_16qam = zeros(length(SNR_Arr), num_bits/4); %scaled noise_16qam initialization
No_16qam = zeros(1,length(SNR_Arr)); %No_16qam initialization for each SNR
Rx_16qam = zeros(length(SNR_Arr), num_bits/4); %Received data after added noise_16qam
%Received bits after comparing to descision regions using decision boundaries
Rx_16qam_bits = zeros(length(SNR_Arr), num_bits);
%an array of vectors that identifies the error_16qam occured in which bits for
%each received stream after adding noise to transmitted stream 
error_16qam = zeros(length(SNR_Arr), num_bits);
BER_16qam = zeros(1,length(SNR_Arr)); %bit error rate for each stream of bits
 
%% 16QAM Mapping
%looping in stream of bits with step = 4 to represent each 4 bits with
%complex numbers to identify constellation points using gray encoding 
for n = 1 : 4 : num_bits
    if [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 0 0 0]
        symbols_16qam(floor(n/4 +1))= -3-3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 0 0 1]
        symbols_16qam(floor(n/4 +1))= -3-1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 0 1 0]
        symbols_16qam(floor(n/4 +1))= -3+3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 0 1 1]
        symbols_16qam(floor(n/4 +1))= -3+1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 1 0 0]
        symbols_16qam(floor(n/4 +1))= -1-3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 1 0 1]
        symbols_16qam(floor(n/4 +1))= -1-1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 1 1 0]
        symbols_16qam(floor(n/4 +1))= -1+3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [0 1 1 1]
        symbols_16qam(floor(n/4 +1))= -1+1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 0 0 0]
        symbols_16qam(floor(n/4 +1))= 3-3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 0 0 1]
        symbols_16qam(floor(n/4 +1))= 3-1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 0 1 0]
        symbols_16qam(floor(n/4 +1))= 3+3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 0 1 1]
        symbols_16qam(floor(n/4 +1))= 3+1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 1 0 0]
        symbols_16qam(floor(n/4 +1))= 1-3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 1 0 1]
        symbols_16qam(floor(n/4 +1))= 1-1j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 1 1 0]
        symbols_16qam(floor(n/4 +1))= 1+3j;
    elseif [data_bits(n) data_bits(n+1) data_bits(n+2) data_bits(n+3)] == [1 1 1 1]
        symbols_16qam(floor(n/4 +1))= 1+1j;
    end
end
 
 
%% Channel
for SNR = SNR_Arr
No_BPSK(SNR) = Esymb_BPSK/(10^(SNR/10)); %No of BPSK for each SNR
No_QPSK(SNR) = (Esymb_QPSK/2)/(10^(SNR/10)); %No of QPSK for each SNR
No_8PSK(SNR) = (Esymb_8PSK/3)/(10^(SNR/10)); %No of 8PSK for each SNR
No_16qam(SNR) = (2.5*Esymb_16qam)/(10^(SNR/10)); %No of 16 QAM for each SNR
 
%making the noise with variance No/2 by multipying it by sqrt(No/2) for
%each modulation scheme
scaled_noise_BPSK(SNR,:) = noise_BPSK*sqrt(No_BPSK(SNR)/2);
scaled_noise_QPSK(SNR,:) = noise_QPSK*sqrt(No_QPSK(SNR)/2);
scaled_noise_8PSK(SNR,:) = noise_8PSK*sqrt(No_8PSK(SNR)/2);
scaled_noise_16qam(SNR,:) = noise_16qam*sqrt(No_16qam(SNR)/2);
 
%Received signal with noise added to it for each modulation scheme
Rx_BPSK(SNR,:) = symbols_BPSK + scaled_noise_BPSK(SNR,:);
Rx_QPSK_GE(SNR,:) = symbols_QPSK_GE + scaled_noise_QPSK(SNR,:);
Rx_QPSK_NGE(SNR,:) = symbols_QPSK_NGE + scaled_noise_QPSK(SNR,:);
Rx_8PSK(SNR,:) = symbols_8PSK + scaled_noise_8PSK(SNR,:);
Rx_16qam(SNR,:) = symbols_16qam + scaled_noise_16qam(SNR,:);
 
%% BPSK Demapping
% transforming received symbols into bits 
for i = 1: num_bits
        if Rx_BPSK(SNR, i) > 0
            Rx_BPSK_bits(SNR, i) = 1 ;
        else
            Rx_BPSK_bits(SNR, i) = 0 ;
        end
end
 
%% QPSK Demapping
%get transmitted bits from received symbols using decision regions of
%constellation points there are 4 decision regions seperated by decision
%boundaries
for i = 1: num_bits/2
    %Grey Encoding Demapping
        if real(Rx_QPSK_GE(SNR, i)) > 0 && imag(Rx_QPSK_GE(SNR, i)) > 0
            Rx_QPSK_GE_bits(SNR, 2*i-1) = 1;
            Rx_QPSK_GE_bits(SNR, 2*i) = 1;
        elseif real(Rx_QPSK_GE(SNR, i)) > 0 && imag(Rx_QPSK_GE(SNR, i)) <= 0
            Rx_QPSK_GE_bits(SNR, 2*i-1) = 1;
            Rx_QPSK_GE_bits(SNR, 2*i) = 0;
        elseif real(Rx_QPSK_GE(SNR, i)) <= 0 && imag(Rx_QPSK_GE(SNR, i)) > 0
            Rx_QPSK_GE_bits(SNR, 2*i-1) = 0;
            Rx_QPSK_GE_bits(SNR, 2*i) = 1;
        elseif real(Rx_QPSK_GE(SNR, i)) <= 0 && imag(Rx_QPSK_GE(SNR, i)) <= 0
            Rx_QPSK_GE_bits(SNR, 2*i-1) = 0;
            Rx_QPSK_GE_bits(SNR, 2*i) = 0;
        end
    
    %Non Grey Encoding Demapping
        if real(Rx_QPSK_NGE(SNR, i)) > 0 && imag(Rx_QPSK_NGE(SNR, i)) > 0
            Rx_QPSK_NGE_bits(SNR, 2*i-1) = 1;
            Rx_QPSK_NGE_bits(SNR, 2*i) = 0;
        elseif real(Rx_QPSK_NGE(SNR, i)) > 0 && imag(Rx_QPSK_NGE(SNR, i)) <= 0
            Rx_QPSK_NGE_bits(SNR, 2*i-1) = 1;
            Rx_QPSK_NGE_bits(SNR, 2*i) = 1;
        elseif real(Rx_QPSK_NGE(SNR, i)) <= 0 && imag(Rx_QPSK_NGE(SNR, i)) > 0
            Rx_QPSK_NGE_bits(SNR, 2*i-1) = 0;
            Rx_QPSK_NGE_bits(SNR, 2*i) = 1;
        elseif real(Rx_QPSK_NGE(SNR, i)) <= 0 && imag(Rx_QPSK_NGE(SNR, i)) <= 0
            Rx_QPSK_NGE_bits(SNR, 2*i-1) = 0;
            Rx_QPSK_NGE_bits(SNR, 2*i) = 0;
        end
end
 
%% 8PSK Demapping
%get transmitted bits from received symbols using decision regions of
%constellation points there are 8 decision regions seperated by decision
%boundaries
for i = 1: num_bits/3
        if  angle(Rx_8PSK(SNR, i)) > -22.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= 22.5*pi/180
            Rx_8PSK_bits(SNR, 3*i-2) = 0;
            Rx_8PSK_bits(SNR, 3*i-1) = 0;
            Rx_8PSK_bits(SNR, 3*i-0) = 0;
        elseif angle(Rx_8PSK(SNR, i)) > 22.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= 67.5*pi/180
            Rx_8PSK_bits(SNR, 3*i-2) = 0;
            Rx_8PSK_bits(SNR, 3*i-1) = 0;
            Rx_8PSK_bits(SNR, 3*i-0) = 1;
        elseif angle(Rx_8PSK(SNR, i)) > 67.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= 112.5*pi/180
            Rx_8PSK_bits(SNR, 3*i-2) = 0;
            Rx_8PSK_bits(SNR, 3*i-1) = 1;
            Rx_8PSK_bits(SNR, 3*i-0) = 1;
        elseif angle(Rx_8PSK(SNR, i)) > 112.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= 157.5*pi/180
            Rx_8PSK_bits(SNR, 3*i-2) = 0;
            Rx_8PSK_bits(SNR, 3*i-1) = 1;
            Rx_8PSK_bits(SNR, 3*i-0) = 0;
        elseif (angle(Rx_8PSK(SNR, i)) > 157.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= 180*pi/180) || (angle(Rx_8PSK(SNR, i)) <= -157.5*pi/180 && angle(Rx_8PSK(SNR, i)) > -180*pi/180)
            Rx_8PSK_bits(SNR, 3*i-2) = 1;
            Rx_8PSK_bits(SNR, 3*i-1) = 1;
            Rx_8PSK_bits(SNR, 3*i-0) = 0;
        elseif angle(Rx_8PSK(SNR, i)) > -157.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= -112.5*pi/180
            Rx_8PSK_bits(SNR, 3*i-2) = 1;
            Rx_8PSK_bits(SNR, 3*i-1) = 1;
            Rx_8PSK_bits(SNR, 3*i-0) = 1;
        elseif angle(Rx_8PSK(SNR, i)) > -112.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= -67.5*pi/180
            Rx_8PSK_bits(SNR, 3*i-2) = 1;
            Rx_8PSK_bits(SNR, 3*i-1) = 0;
            Rx_8PSK_bits(SNR, 3*i-0) = 1;
        elseif angle(Rx_8PSK(SNR, i)) > -67.5*pi/180 && angle(Rx_8PSK(SNR, i)) <= -22.5*pi/180
            Rx_8PSK_bits(SNR, 3*i-2) = 1;
            Rx_8PSK_bits(SNR, 3*i-1) = 0;
            Rx_8PSK_bits(SNR, 3*i-0) = 0;
        end
end
 
%% 16QAM Demapping
%get transmitted bits from received symbols using decision regions of
%constellation points there are 16 decision regions seperated by decision
%boundaries
for i = 1: num_bits/4
    if  real(Rx_16qam(SNR, i)) <= -2 && imag(Rx_16qam(SNR, i)) <= -2
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) <= -2 && imag(Rx_16qam(SNR, i)) < 0 && imag(Rx_16qam(SNR, i)) > -2
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    elseif real(Rx_16qam(SNR, i)) <= -2 && imag(Rx_16qam(SNR, i)) >= 2
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) <= -2 && imag(Rx_16qam(SNR, i)) < 2 && imag(Rx_16qam(SNR, i)) >= 0
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    elseif real(Rx_16qam(SNR, i)) < 0 && real(Rx_16qam(SNR, i)) > -2 && imag(Rx_16qam(SNR, i)) <= -2
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) < 0 && real(Rx_16qam(SNR, i)) > -2 && imag(Rx_16qam(SNR, i)) > -2 && imag(Rx_16qam(SNR, i)) < 0
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    elseif real(Rx_16qam(SNR, i)) < 0 && real(Rx_16qam(SNR, i)) > -2 && imag(Rx_16qam(SNR, i)) >= 2
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) < 0 && real(Rx_16qam(SNR, i)) > -2 && imag(Rx_16qam(SNR, i)) >= 0 && imag(Rx_16qam(SNR, i)) < 2
        Rx_16qam_bits(SNR, 4*i-3) = 0;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    elseif  real(Rx_16qam(SNR, i)) >= 2 && imag(Rx_16qam(SNR, i)) <= -2
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) >= 2 && imag(Rx_16qam(SNR, i)) < 0 && imag(Rx_16qam(SNR, i)) > -2
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    elseif real(Rx_16qam(SNR, i)) >= 2 && imag(Rx_16qam(SNR, i)) >= 2
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) >= 2 && imag(Rx_16qam(SNR, i)) < 2 && imag(Rx_16qam(SNR, i)) >= 0
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 0;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    elseif real(Rx_16qam(SNR, i)) >= 0 && real(Rx_16qam(SNR, i)) < 2 && imag(Rx_16qam(SNR, i)) <= -2
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) >= 0 && real(Rx_16qam(SNR, i)) < 2 && imag(Rx_16qam(SNR, i)) > -2 && imag(Rx_16qam(SNR, i)) < 0
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 0;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    elseif real(Rx_16qam(SNR, i)) >= 0 && real(Rx_16qam(SNR, i)) < 2 && imag(Rx_16qam(SNR, i)) >= 2
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 0;
    elseif real(Rx_16qam(SNR, i)) >= 0 && real(Rx_16qam(SNR, i)) < 2 && imag(Rx_16qam(SNR, i)) >= 0 && imag(Rx_16qam(SNR, i)) < 2
        Rx_16qam_bits(SNR, 4*i-3) = 1;
        Rx_16qam_bits(SNR, 4*i-2) = 1;
        Rx_16qam_bits(SNR, 4*i-1) = 1;
        Rx_16qam_bits(SNR, 4*i-0) = 1;
    end
end
 
%% Error Calculation 
%recognize the bits that differed between transmitted and received streams 
error_BPSK(SNR,:) = Rx_BPSK_bits(SNR,:) == data_bits;
error_QPSK_GE(SNR,:) = Rx_QPSK_GE_bits(SNR,:) == data_bits;
error_QPSK_NGE(SNR,:) = Rx_QPSK_NGE_bits(SNR,:) == data_bits;
error_8PSK(SNR,:) = Rx_8PSK_bits(SNR,:) == data_bits;
error_16qam(SNR,:) = Rx_16qam_bits(SNR,:) == data_bits;
 
%counting the number of flipped bits and get the ratio
%between number of bits that was flipped and the total number of bits in
%each stream for each SNR to get bit error rate for each SNR value for each
%modulation scheme
BER_BPSK(SNR) = num_bits - nnz(error_BPSK(SNR,:));
BER_BPSK(SNR) = BER_BPSK(SNR)/num_bits;
 
BER_QPSK_GE(SNR) = num_bits - nnz(error_QPSK_GE(SNR,:));
BER_QPSK_GE(SNR) = BER_QPSK_GE(SNR)/num_bits;
 
BER_QPSK_NGE(SNR) = num_bits - nnz(error_QPSK_NGE(SNR,:));
BER_QPSK_NGE(SNR) = BER_QPSK_NGE(SNR)/num_bits;
 
BER_8PSK(SNR) = num_bits - nnz(error_8PSK(SNR,:));
BER_8PSK(SNR) = BER_8PSK(SNR)/num_bits;
 
BER_16qam(SNR) = num_bits - nnz(error_16qam(SNR,:));
BER_16qam(SNR) = BER_16qam(SNR)/num_bits;
end
 
%% Theoritical Calculations
%theoretical value of bit error rate for each modulation scheme
BER_BPSK_Theoritical = 0.5*erfc(sqrt(Esymb_BPSK./No_BPSK));
BER_QPSK_Theoritical = 0.5*erfc(sqrt((Esymb_QPSK/2)./No_QPSK));
BER_8PSK_Theoritical = erfc(sqrt(Esymb_8PSK./No_8PSK)*sin(pi/8))/3;
BER_16qam_Theoritical = 1.5*erfc(sqrt(Esymb_16qam./No_16qam))/4;
 
%plotting bit error rate with log scale versus SNR for all modulation schemes and
%drawing them on the same graph along with the theoretical BER of all
%modulation schemes
figure('Name', 'BER Vs SNR');
semilogy(SNR_Arr, BER_BPSK_Theoritical, '--r', SNR_Arr, BER_BPSK, 'c', 'LineWidth', 1.5);
hold on
semilogy(SNR_Arr, BER_QPSK_Theoritical, '--r', SNR_Arr, BER_QPSK_GE, 'g' ,'LineWidth', 1.5);
hold on
semilogy(SNR_Arr, BER_8PSK_Theoritical, '--m', SNR_Arr, BER_8PSK, 'b', 'LineWidth', 1.5);
hold on
semilogy(SNR_Arr, BER_16qam_Theoritical, '--y', SNR_Arr, BER_16qam, 'k', 'LineWidth', 1.5);
title ("BER Vs SNR");
xlabel("Eb/No (db)");
ylabel("BER");
xlim([1,10]);
legend ('BER BPSK Theoretical', 'BER BPSK Practical', 'BER QPSK Theoretical', 'BER QPSK Practical', 'BER 8PSK Theoretical', 'BER 8PSK Practical', 'BER 16qam Theoretical', 'BER 16qam Practical');
grid on;
 
%% Comparison between gray and non gray encoding in QPSK modulation
figure('Name', 'BER Vs SNR');
semilogy(SNR_Arr, BER_QPSK_GE, 'r', SNR_Arr, BER_QPSK_NGE, 'k', 'LineWidth', 1.5);
title ("BER Vs SNR");
xlabel("Eb/N0 (db)");
ylabel("BER");
xlim([1,10]);
ylim([10^-3,1]);
legend ('BER QPSK Gray Encoding', 'BER QPSK Not Gray Encoding')
grid on;
