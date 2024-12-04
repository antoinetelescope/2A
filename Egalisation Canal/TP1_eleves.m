






% NOTES : 2.1)  lorsqu'on met le dsp de z + n ds la boucle en ii, on voit
% que plus le snr est petit (DONC TEB GRAND), plusla dsp se rapproche de celle d'un bruit blanc (constante)
% Quand le snr est plus grand (dernières itérations), la dsp ressemble
% davantage à celle de h (y = hx et dsp de x = cste).





% Script for computing the BER for BPSK/QPSK modulation in ISI Channels
% 
close all;
clear all;

%% Simulation parameters
% On décrit ci-après les paramètres généraux de la simulation

%Frame length
M=4; %2:BPSK, 4: QPSK
N  = 10000; % Number of transmitted bits or symbols
Es_N0_dB = [0:3:30]; % Eb/N0 values
%Multipath channel parameters
hc=[1 0.8*exp(1i*pi/3) 0.3*exp(1i*pi/6) ];%0.1*exp(1i*pi/12)];%ISI channel
% hc=[0.04, -0.05, 0.07, -0.21, -0.5, 0.72, 0.36, 0, 0.21, 0.03, 0.07];
%hc=[0.407, 0.815, 0.407]
%hc = [0.227, 0.46, 0.688, 0.460, 0.227];
%a=0.2;
%hc=[1 -a];
Lc=length(hc);%Channel length
ChannelDelay=0; %delay is equal to number of non causal taps
%Preallocations
nErr_zfinf=zeros(1,length(Es_N0_dB));
for ii = 1:length(Es_N0_dB)

   % BPSK symbol generations
%    bits = rand(1,N)>0.5; % generating 0,1 with equal probability
%    s = 1-2*bits; % BPSK modulation following: {0 -> +1; 1 -> -1} 
   
    % QPSK symbol generations
   bits = rand(2,N)>0.5; % generating 0,1 with equal probability
   s = 1/sqrt(2)*((1-2*bits(1,:))+1j*(1-2*bits(2,:))); % QPSK modulation following the BPSK rule for each quadatrure component: {0 -> +1; 1 -> -1} 
   sigs2=var(s);
   
   % Channel convolution: equivalent symbol based representation
   z = conv(hc,s);
   
   %Generating noise
   sig2b=10^(-Es_N0_dB(ii)/10);
   %n = sqrt(sig2b)*randn(1,N+Lc-1); % white gaussian noise, BPSK Case
    n = sqrt(sig2b/2)*randn(1,N+Lc-1)+1j*sqrt(sig2b/2)*randn(1,N+Lc-1); % white gaussian noise, QPSK case
   
   % Adding Noise

   y = z + n; % additive white gaussian noise
   u = z;


   %% zero forcing equalization
   % We now study ZF equalization
   
   %Unconstrained ZF equalization, only if stable inverse filtering
   
   
   %%
   % 
   %  The unconstrained ZF equalizer, when existing is given by 
   % 
   % $w_{,\infty,zf}=\frac{1}{h(z)}$
   % 
   %%
   % 
   
    s_zf=filter(1,hc,y);%if stable causal filter is existing
    bhat_zf = zeros(2,length(bits));
    bhat_zf(1,:)= real(s_zf(1:N)) < 0;
    bhat_zf(2,:)= imag(s_zf(1:N)) < 0;
    nErr_zfinfdirectimp(1,ii) = size(find([bits(:)-bhat_zf(:)]),1);

    %Otherwise, to handle the non causal case
    Nzf=200;
    [r, p, k]=residuez(1, hc);
    [w_zfinf]=ComputeRI( Nzf, r, p, k );
    s_zf=conv(w_zfinf,y);
    bhat_zf = zeros(2,length(bits));
    bhat_zf(1,:)= real(s_zf(Nzf:N+Nzf-1)) < 0;
    bhat_zf(2,:)= imag(s_zf(Nzf:N+Nzf-1)) < 0;
    nErr_zfinf(1,ii) = size(find([bits(:)- bhat_zf(:)]),1);
    
    %%MMSE

    deltac= zeros(1,2*Lc-1);
    deltac(Lc) = 1 ;
    Nmmse = 200;% non causal
    [r,p,k]=residuez(fliplr(conj(hc)),(conv(hc,fliplr(conj(hc)))+(sig2b/sigs2)*deltac));
    [w_mmseinf]=ComputeRI(Nmmse,r,p,k);
    s_mmse=conv(w_mmseinf,y);

    bhat_mmse = zeros(2,length(bits));
    bhat_mmse(1,:)= real(s_mmse(Nmmse:N+Nmmse-1)) < 0;
    bhat_mmse(2,:)= imag(s_mmse(Nmmse:N+Nmmse-1)) < 0;
    nErr_mmseinfdirectimp(1,ii) = size(find([bits(:)- bhat_mmse(:)]),1);

    nErr_mmseinf(1,ii) = size(find([bits(:)- bhat_mmse(:)]),1);




   % DSP3 = pwelch(s_mmse,[],[],[],'centered');
   % figure('Name','DSP')
   % plot(10*log(DSP3))
   % xlabel('Fréquence (Hz)')
   % ylabel('DSP Sortie MMSE')


   %% FIR ZF

    Nw=12;
    d = 5;          %distance
    H = toeplitz([hc(1) zeros(1,Nw-1)]',[hc, zeros(1,Nw-1)]);
    Ry = (conj(H)*H.');
    p = zeros(Nw+Lc-1,1);

    P= (H.'*inv((Ry))*conj(H));
    [alpha,dopt]=max(diag(abs(P)));
    %p(d+1)=1
    p(dopt)=1;
    Gamma = conj(H)*p;
    w_zf_fir = (inv(Ry)*Gamma).';

    sig_e_opt = sigs2-conj(w_zf_fir)*Gamma;
    bias = 1-sig_e_opt/sigs2;
    shat = conv(w_zf_fir,y);
    shat = shat(dopt:end);

    bHat = zeros(2,length(bits));
    bHat(1,:)=real(shat(1:N))<0;
    bHat(2,:)=imag(shat(1:N))<0;

    nErr_zf_ls(1,ii) = size(find([bits(:)- bHat(:)]),1);

    %% FIR MMSE

    Nw = 12;
    d=5;
    H = toeplitz([hc(1) zeros(1,Nw-1)]',[hc, zeros(1,Nw-1)]);
    Ry = sigs2*(conj(H)*H.')+sig2b*eye(Nw);
    p = zeros(Nw+Lc-1,1);

    P = 1/sigs2*(H.'*inv((Ry/sigs2))*conj(H));
    [alpha,dopt]=max(diag(abs(P)));

    %p(d+1)=1
    p(dopt)=1;

    Gamma = sigs2*conj(H)*p;
    w_mmse_fir = (inv(Ry)*Gamma).';
    
    sig_e_opt = sigs2 -conj(w_mmse_fir)*Gamma;
    bias = 1-sig_e_opt/sigs2;
    shat = conv(w_mmse_fir,y);
    shat = shat(dopt:end);

    bHat = zeros(2,length(bits));
    bHat(1,:)=real(shat(1:N))<0;
    bHat(2,:)=imag(shat(1:N))<0;

    nErr_mmse_ls(1,ii) = size(find([bits(:)- bHat(:)]),1);

    %% ML

    % s_ml = mlseeq(y,hc,const,tblen,'rst',nsamp,[],[]) ;

end

simBer_mmseinf = nErr_mmseinf/N/log2(M);
simBer_mmseinfdirectimp = nErr_mmseinfdirectimp/N/log2(M);

simBer_zfinfdirectimp = nErr_zfinfdirectimp/N/log2(M); 
simBer_zfinf = nErr_zfinf/N/log2(M); 

simBer_zf_ls= nErr_zf_ls/N/log2(M);

simBer_mmse_ls= nErr_mmse_ls/N/log2(M);


% plot

figure


semilogy(Es_N0_dB,simBer_zfinfdirectimp(1,:),'p-','Linewidth',2);
% hold on
% semilogy(Es_N0_dB,simBer_zfinf(1,:),'bs-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_mmseinfdirectimp(1,:),'p-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_zf_ls(1,:),'p-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_mmse_ls(1,:),'p-','Linewidth',2);
% hold on
% semilogy(Es_N0_dB,simBer_mmseinf(1,:),'rs-','Linewidth',2);
axis([0 40 10^-6 0.5])
grid on
legend('sim-zf-inf/direct','sim-mmse-inf/direct','sim-zf-rif','sim-mmse-rif');
xlabel('E_s/N_0, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK in ISI with ZF and MMSE equalizers')

figure
title('Impulse response')
stem(real(w_mmseinf))
hold on
stem(real(w_zfinf),'r-')
ylabel('Amplitude');
xlabel('time index')


% DSP = pwelch(y,[],[],[],'centered');
% figure('Name','DSP')
% plot(10*log(DSP))
% xlabel('Fréquence (Hz)')
% ylabel('DSP avec bruit')
% 
% DSP1 = pwelch(u,[],[],[],'centered');
% figure('Name','DSP')
% plot(10*log(DSP1))
% xlabel('Fréquence (Hz)')
% ylabel('DSP sans bruit')


% DSP2 = pwelch(s_zf,[],[],[],'centered');
% figure('Name','DSP')
% plot(10*log(DSP2))
% xlabel('Fréquence (Hz)')
% ylabel('DSP Sortie ZF')







