
%% Initialize
SNR=1:1:40;                 %Range of SNR in dB
N=2^14;                  %Number of bits for train
target_BER = 3e-3;       %Define the thresholds
max_iter=10;                %Define the maximum iteration
Mo1=2;                       %BPSK
Mo2=4;                       %QPSK
Mo3=16;                      %16QAM
thresholds1=20;
thresholds2=23;
thresholds3=29;
Data_BPSK=[];
Data_QPSK=[];
Data_16QAM=[];
Output_Data_BPSK_BER=[];
Output_Data_BPSK_SNR=[];
Output_Data_BPSK=[];
Output_Data_QPSK_BER=[];
Output_Data_QPSK_SNR=[];
Output_Data_QPSK=[];
Output_Data_16QAM_BER=[];
Output_Data_16QAM_SNR=[];
Output_Data_16QAM=[];
Output_Data_BER_SNR=[];
Output_Data_Class=[];
KNN_BPSK_counter=[];
KNN_QPSK_counter=[];
KNN_16QAM_counter=[];

%% Generate Train Data
for iter=1:1:max_iter
    
    if iter<=7
        R=raylrnd(0.5,1,N);         %Produce the Rayleigh signal
    else
        R=raylrnd(0.8,1,N);         %Produce another Rayleigh signal
    end
    
    x1=randi([0,Mo1-1],1,N);     %Produce the random signal
    x2=randi([0,Mo2-1],1,N);
    x3=randi([0,Mo3-1],1,N);
    h1=pskmod(x1,Mo1);            %BPSK Modulation
    h2=pskmod(x2,Mo2);            %QPSK Modulation
    h3=qammod(x3,Mo3);            %16QAM Modulation
    H1=h1.*R;                   %BPSK with Rayleigh Channel
    H2=h2.*R;                   %QPSK with Rayleigh Channel
    H3=h3.*R;                   %16QAM with Rayleigh Channel

    for i=1:length(SNR)

        y_RE_n1_Train=R.\awgn(H1,SNR(i),'measured');
        y_RE_1_Train=pskdemod(y_RE_n1_Train,Mo1);     
        [bit_RE1,ratio1]=biterr(x1,y_RE_1_Train);
        BPSK_Ray_Train(i)=ratio1;

        y_RE_n2_Train=R.\awgn(H2,SNR(i),'measured');
        y_RE_2_Train=pskdemod(y_RE_n2_Train,Mo2);     
        [bit_RE2,ratio2]=biterr(x2,y_RE_2_Train);
        QPSK_Ray_Train(i)=ratio2;

        y_RE_n3_Train=R.\awgn(H3,SNR(i),'measured');
        y_RE_3_Train=qamdemod(y_RE_n3_Train,Mo3);     
        [bit_RE3,ratio3]=biterr(x3,y_RE_3_Train);
        QAM_Ray_Train(i)=ratio3;

    end
    Data_BPSK=[BPSK_Ray_Train;Data_BPSK];
    Data_QPSK=[QPSK_Ray_Train;Data_QPSK];
    Data_16QAM=[QAM_Ray_Train;Data_16QAM];

end

%% Plot
figure(1)
semilogy(SNR,Data_BPSK,'ro');%BPSK
hold on;
grid on; 
semilogy(SNR,Data_QPSK,'go');%QPSK
semilogy(SNR,Data_16QAM,'bo');%16QAM
title('Distribution of Training Data');
xlabel('SNR（dB）');ylabel('BER');
line([0 40],[target_BER target_BER],'Color','red','LineStyle','--')
hold off;

%% Reform and process data

for i=1:40
    Output_Data_BPSK_BER=[Output_Data_BPSK_BER;Data_BPSK(:,i)];
    Output_Data_BPSK_SNR=[Output_Data_BPSK_SNR;i.*ones(10,1)];
    Output_Data_QPSK_BER=[Output_Data_QPSK_BER;Data_QPSK(:,i)];
    Output_Data_QPSK_SNR=[Output_Data_QPSK_SNR;i.*ones(10,1)];
    Output_Data_16QAM_BER=[Output_Data_16QAM_BER;Data_16QAM(:,i)];
    Output_Data_16QAM_SNR=[Output_Data_16QAM_SNR;i.*ones(10,1)];
end
%Output it in BER,SNR ~ CLASS style
Output_Data_BPSK=[Output_Data_BPSK_BER,Output_Data_BPSK_SNR];
Output_Data_QPSK=[Output_Data_QPSK_BER,Output_Data_QPSK_SNR];
Output_Data_16QAM=[Output_Data_16QAM_BER,Output_Data_16QAM_SNR];
Output_Data_BER_SNR=[Output_Data_BPSK;Output_Data_QPSK;Output_Data_16QAM];
Output_Data_Class=[repmat({'BPSK'},400,1);repmat({'QPSK'},400,1);repmat({'16QAM'},400,1)];

index=find(Output_Data_BER_SNR(:,1)<=target_BER);
Require_Train_Data=Output_Data_BER_SNR(index,:);
Require_Train_Class=Output_Data_Class(index);

%% Train Data and make predictions on different SNR
KNNC = fitcknn(Require_Train_Data,Require_Train_Class,'NumNeighbors',20,'Standardize',1)
New_Sample=[]

for Pred_SNR=1:1:40
    %Make predictions raw data set

    New_Sample = [New_Sample;target_BER Pred_SNR];
end
%Get predictions
[label,score,cost] = predict(KNNC,New_Sample);

%% Generate Data
Practical_SNR=thresholds1:1:40;                 %Range of SNR in dB
N=2^20;                  %Number of bits to simulate
x1=randi([0,Mo1-1],1,N);     %Produce the random signal
x2=randi([0,Mo2-1],1,N);
x3=randi([0,Mo3-1],1,N);
R=raylrnd(0.5,1,N);         %Produce the Rayleigh signal

%% Modulation and transmit over a Rayleigh-distributed narrowband fading channel.
h1=pskmod(x1,Mo1);            %BPSK Modulation
h2=pskmod(x2,Mo2);            %QPSK Modulation
h3=qammod(x3,Mo3);            %16QAM Modulation
H1=h1.*R;                   %BPSK with Rayleigh Channel
H2=h2.*R;                   %QPSK with Rayleigh Channel
H3=h3.*R;                   %16QAM with Rayleigh Channel

%% KNN-AWGN and Demodulation
for Practical_SNR=thresholds1:1:40
    mod=label{Practical_SNR}
    if strcmp(mod,'BPSK')
        y_RE_n1=R.\awgn(H1,Practical_SNR,'measured');
        y_RE_1=pskdemod(y_RE_n1,Mo1);     
        [bit_RE1,ratio]=biterr(x1,y_RE_1);
        KNN_BER(Practical_SNR-(thresholds1-1))=ratio;
        KNN_BPSK_counter=[KNN_BPSK_counter;1];
        continue
    elseif strcmp(mod,'QPSK')
        y_RE_n2=R.\awgn(H2,Practical_SNR,'measured');
        y_RE_2=pskdemod(y_RE_n2,Mo2);     
        [bit_RE2,ratio]=biterr(x2,y_RE_2);
        if ratio > target_BER
            y_RE_n1=R.\awgn(H1,Practical_SNR,'measured');
            y_RE_1=pskdemod(y_RE_n1,Mo1);     
            [bit_RE1,ratio]=biterr(x1,y_RE_1);
            KNN_BER(Practical_SNR-(thresholds1-1))=ratio;
            KNN_BPSK_counter=[KNN_BPSK_counter;1];
            continue
        end
        KNN_BER(Practical_SNR-(thresholds1-1))=ratio;
        KNN_QPSK_counter=[KNN_QPSK_counter;1];
        continue
    elseif strcmp(mod,'16QAM')
        y_RE_n3=R.\awgn(H3,Practical_SNR,'measured');
        y_RE_3=qamdemod(y_RE_n3,Mo3);     
        [bit_RE3,ratio]=biterr(x3,y_RE_3);
        if ratio > target_BER
            y_RE_n2=R.\awgn(H2,Practical_SNR,'measured');
            y_RE_2=pskdemod(y_RE_n2,Mo2);     
            [bit_RE2,ratio]=biterr(x2,y_RE_2);
            KNN_BER(Practical_SNR-(thresholds1-1))=ratio;
            KNN_QPSK_counter=[KNN_QPSK_counter;1];
            continue
        end
        KNN_BER(Practical_SNR-(thresholds1-1))=ratio;
        KNN_16QAM_counter=[KNN_16QAM_counter;1];
        continue
    end
end

%% KNN Throughput
NO_trans = 0*ones(thresholds1-1,1);
KNN_BPSK_trans = 2*KNN_BPSK_counter;
KNN_QPSK_trans = 4*KNN_QPSK_counter;
KNN_16QAM_trans = 16*KNN_16QAM_counter;
KNN_Throughput = [NO_trans;KNN_BPSK_trans;KNN_QPSK_trans;KNN_16QAM_trans];

%% AWGN and Demodulation for each modulation
for Practical_SNR=thresholds1:1:40

    y_RE_n1=R.\awgn(H1,Practical_SNR,'measured');
    y_RE_1=pskdemod(y_RE_n1,Mo1);     
    [bit_RE1,ratio1]=biterr(x1,y_RE_1);
    BPSK_Ray_BER(Practical_SNR-(thresholds1-1))=ratio1;
    
    y_RE_n2=R.\awgn(H2,Practical_SNR,'measured');
    y_RE_2=pskdemod(y_RE_n2,Mo2);     
    [bit_RE2,ratio2]=biterr(x2,y_RE_2);
    QPSK_Ray_BER(Practical_SNR-(thresholds1-1))=ratio2;
    
    y_RE_n3=R.\awgn(H3,Practical_SNR,'measured');
    y_RE_3=qamdemod(y_RE_n3,Mo3);     
    [bit_RE3,ratio3]=biterr(x3,y_RE_3);
    QAM_Ray_BER(Practical_SNR-(thresholds1-1))=ratio3;
    
end

%% Adaptive Demodulation
for Practical_SNR=thresholds1:1:40
    if Practical_SNR>=thresholds1 && Practical_SNR<thresholds2
        y_RE_n1=R.\awgn(H1,Practical_SNR,'measured');
        y_RE_1=pskdemod(y_RE_n1,Mo1);     
        [bit_RE1,ratio]=biterr(x1,y_RE_1);
        Adaptive_Modulation(Practical_SNR-(thresholds1-1))=ratio;
    elseif Practical_SNR>=thresholds2 && Practical_SNR<thresholds3
        y_RE_n2=R.\awgn(H2,Practical_SNR,'measured');
        y_RE_2=pskdemod(y_RE_n2,Mo2);     
        [bit_RE2,ratio]=biterr(x2,y_RE_2);
        Adaptive_Modulation(Practical_SNR-(thresholds1-1))=ratio;
    elseif Practical_SNR>=thresholds3
        y_RE_n3=R.\awgn(H3,Practical_SNR,'measured');
        y_RE_3=qamdemod(y_RE_n3,Mo3);     
        [bit_RE3,ratio]=biterr(x3,y_RE_3);
        Adaptive_Modulation(Practical_SNR-(thresholds1-1))=ratio;
    end
end

%% Adaptive Modulation Throughput

BPSK_trans = 2*ones(thresholds2-thresholds1,1);
QPSK_trans = 4*ones(thresholds3-thresholds2,1);
QAM_trans = 16*ones(40-thresholds3+1,1);
Adaptive_Modulation_Throughput = [NO_trans;BPSK_trans;QPSK_trans;QAM_trans];

%% Plot figure
figure(2)
Practical_SNR=thresholds1:1:40; 
semilogy(Practical_SNR,BPSK_Ray_BER,':rx');
hold on;
semilogy(Practical_SNR,QPSK_Ray_BER,':gx');
semilogy(Practical_SNR,QAM_Ray_BER,':bx');
semilogy(Practical_SNR,KNN_BER,':ko'); 
semilogy(Practical_SNR,Adaptive_Modulation,':m*');

grid on;
axis([15 40 10^-5 1.2]);
line([0 40],[target_BER target_BER],'Color','red','LineStyle','--')
legend('BPSK','QPSK','16QAM','KNN aided adaptive','Original Adaptive');
title('KNN aided adaptive system (BER vs SNR)');
xlabel('SNR（dB）');ylabel('BER');

%% Plot
figure(3)
axis([0 40 0 20]);
plot(SNR,Adaptive_Modulation_Throughput,'Linewidth',2);
hold on;
plot(SNR,KNN_Throughput,'Linewidth',2,'LineStyle','--');

legend({'Adaptive Modulation','KNN Throughput'},'Location','northwest');
title('Throughput vs SNR');
xlabel('SNR（dB）');ylabel('Throughput');
hold off;





