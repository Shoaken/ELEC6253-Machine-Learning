
%% initialize Data
N=2^20;                  %Number of bits to simulate
target_BER = 3e-3;       %Define the thresholds
Mo1=2;                       %BPSK
Mo2=4;                       %QPSK
Mo3=16;                      %16QAM
x1=randi([0,Mo1-1],1,N);     %Produce the random signal
x2=randi([0,Mo2-1],1,N);
x3=randi([0,Mo3-1],1,N);
R=raylrnd(0.5,1,N);         %Produce the Rayleigh signal
Adaptive_Modulation_THROUGHPUT=[];
thresholds1=20;
thresholds2=23;
thresholds3=29;

%% Modulation and transmit over a Rayleigh-distributed narrowband fading channel.
h1=pskmod(x1,Mo1);            %BPSK Modulation
h2=pskmod(x2,Mo2);            %QPSK Modulation
h3=qammod(x3,Mo3);            %16QAM Modulation
H1=h1.*R;                   %BPSK with Rayleigh Channel
H2=h2.*R;                   %QPSK with Rayleigh Channel
H3=h3.*R;                   %16QAM with Rayleigh Channel

%% AWGN and Demodulation for each modulation
for SNR=1:1:40

    y_RE_n1=R.\awgn(H1,SNR,'measured');
    y_RE_1=pskdemod(y_RE_n1,Mo1);     
    [bit_RE1,ratio1]=biterr(x1,y_RE_1);
    BPSK_Ray(SNR)=ratio1;
    
    y_RE_n2=R.\awgn(H2,SNR,'measured');
    y_RE_2=pskdemod(y_RE_n2,Mo2);     
    [bit_RE2,ratio2]=biterr(x2,y_RE_2);
    QPSK_Ray(SNR)=ratio2;
    
    y_RE_n3=R.\awgn(H3,SNR,'measured');
    y_RE_3=qamdemod(y_RE_n3,Mo3);     
    [bit_RE3,ratio3]=biterr(x3,y_RE_3);
    QAM_Ray(SNR)=ratio3;
    
end

%% Adaptive Demodulation
for SNR=1:1:40
    if SNR <thresholds1                       %Choose Modulation method
       continue;                     %Don't transmit data
    end
    if SNR>=thresholds1 && SNR<thresholds2
        y_RE_n1=R.\awgn(H1,SNR,'measured');
        y_RE_1=pskdemod(y_RE_n1,Mo1);     
        [bit_RE1,ratio]=biterr(x1,y_RE_1);
        Adaptive_Modulation(SNR)=ratio;
    elseif SNR>=thresholds2 && SNR<thresholds3
        y_RE_n2=R.\awgn(H2,SNR,'measured');
        y_RE_2=pskdemod(y_RE_n2,Mo2);     
        [bit_RE2,ratio]=biterr(x2,y_RE_2);
        Adaptive_Modulation(SNR)=ratio;
    elseif SNR>=thresholds3
        y_RE_n3=R.\awgn(H3,SNR,'measured');
        y_RE_3=qamdemod(y_RE_n3,Mo3);     
        [bit_RE3,ratio]=biterr(x3,y_RE_3);
        Adaptive_Modulation(SNR)=ratio;
    end
end

%% Plot figure
figure(1)
SNR=1:1:40;                 %Range of SNR in dB
axis([1 40 10^-5 1]);
semilogy(SNR,BPSK_Ray,':rx');
hold on;
semilogy(SNR,QPSK_Ray,':gx');
semilogy(SNR,QAM_Ray,':bx');
grid on;
line([0 40],[target_BER target_BER],'Color','red','LineStyle','--')
semilogy(thresholds1,BPSK_Ray(thresholds1),'ko');
semilogy(thresholds2,QPSK_Ray(thresholds2),'ko');
semilogy(thresholds3,QAM_Ray(thresholds3),'ko');
legend('BPSK Ray Equalize','QPSK Ray Equalize','16QAM Ray Equalize');
title('the BER curves for each modulation technique');
xlabel('SNR（dB）');ylabel('BER');
hold off;

figure(2)
SNR=1:1:40;                 %Range of SNR in dB
axis([1 40 10^-5 1]);
semilogy(SNR,BPSK_Ray,':rx');
hold on;
semilogy(SNR,QPSK_Ray,':gx');
semilogy(SNR,QAM_Ray,':bx');
grid on;
line([0 40],[target_BER target_BER],'Color','red','LineStyle','--')
semilogy(SNR,Adaptive_Modulation,':ko');
legend('BPSK Ray','QPSK Ray','16QAM Ray','Adaptive Modulation');
title('adaptive modulation system-BER vs SNR');
xlabel('SNR（dB）');ylabel('BER');
hold off;

%% Throughput
NO_trans = 0*ones(thresholds1-1,1);
BPSK_trans = 2*ones(thresholds2-thresholds1,1);
QPSK_trans = 4*ones(thresholds3-thresholds2,1);
QAM_trans = 16*ones(40-thresholds3+1,1);
Adaptive_Modulation_THROUGHPUT = [NO_trans;BPSK_trans;QPSK_trans;QAM_trans];

Fixed_BPSK_trans = 2*ones(40-thresholds1+1,1);
Fixed_QPSK_trans = 4*ones(40-thresholds2+1,1);
Fixed_QPSK_NO_trans = 0*ones(thresholds2-1,1);
Fixed_QAM_trans = 16*ones(40-thresholds3+1,1);
Fixed_QAM_NO_trans = 0*ones(thresholds3-1,1);

Fixed_BPSK_THROUGHPUT = [NO_trans;Fixed_BPSK_trans];
Fixed_QPSK_THROUGHPUT = [Fixed_QPSK_NO_trans;Fixed_QPSK_trans];
Fixed_QAM_THROUGHPUT = [Fixed_QAM_NO_trans;Fixed_QAM_trans];

%% Plot
figure(3)
axis([0 40 0 20]);
plot(SNR,Adaptive_Modulation_THROUGHPUT,'Linewidth',2);
hold on;
plot(SNR,Fixed_BPSK_THROUGHPUT,'--',SNR,Fixed_QPSK_THROUGHPUT,'--',SNR,Fixed_QAM_THROUGHPUT,'--');
legend({'Adaptive Modulation','Fixed BPSK','Fixed QPSK','Fixed 16QAM'},'Location','northwest');
title('Throughput vs SNR');
xlabel('SNR（dB）');ylabel('Throughput');
hold off;









