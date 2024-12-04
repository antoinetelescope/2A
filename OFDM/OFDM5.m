close all;
clear all;

% Variables d'initialisation

N = 16;     % Nombre de porteuses
Nbits = 1000;   % Nombre de bits total
Nactives = 16;  % Nombre de porteuses actives
lignes_garde = 2;   % Nombre de lignes de garde
EBN0dB = 0:0.01:6;  % Valeurs de Eb/N0 en dB
EBN0 = 10.^(EBN0dB/10);    % Valeurs de Eb/N0
M = 4;           % Modulation QPSK

% Génération des bits

bits = randi([0, 1], 1, N*Nbits*2);

% Mapping QPSK
Reel = 1-2*bits(1:2:end);          % Partie réelle
Imaginaire = 1-2*bits(2:2:end);    % Partie imaginaire
Mapping = Reel + 1i*Imaginaire; % Symboles normalisés

MappingMatrice = reshape(Mapping,N,Nbits);

% Canal AWGN

Xe = ifft(MappingMatrice,N);
garde = Xe(N-lignes_garde+1:end,:); 
Xe_garde = [garde;Xe];  % On introduit ici un préfixe cyclique "au dessus" de Xe

ofdm_lineaire = reshape(Xe_garde, 1, []);

Px = mean(abs(ofdm_lineaire).^2);

for i=1:length(EBN0dB)      % 61 : Nombre de points de EBN0dB

    sigmacarre = Px/(2*log2(M)*EBN0(i)); % Puissance du bruit
    noise = sqrt(sigmacarre)*(randn(1, length(ofdm_lineaire)) + 1i*randn(1, length(ofdm_lineaire))); % Bruit complexe
    signalrecu = ofdm_lineaire + noise;
   
    % Démodulation

    Y_reshape = reshape(signalrecu,size(Xe_garde));      % On reshape la matrice en ligne
    Xe_sansprefixe = Y_reshape(lignes_garde+1:N+lignes_garde,:);      % On supprime le préfixe cyclique
    Y_recu = fft(Xe_sansprefixe,N);
    Yligne = reshape(Y_recu,1,[]);

    % Démapping
    
    % Simulation d'une transmission parfaite (sans bruit)
    received_reel = real(Yligne) >= 0;      % Partie réelle reçue
    received_imaginaire = imag(Yligne) >= 0; % Partie imaginaire reçue

    % Reconstruction des bits
    received_bits = zeros(1, length(bits));
    received_bits(1:2:end) = ~received_reel;        % Reconstruction des bits réels
    received_bits(2:2:end) = ~received_imaginaire;  % Reconstruction des bits imaginaires

    TEB(i) = mean(bits ~= received_bits);
end

% Tracé du TEB

TEB_theorique = 2*qfunc(sqrt(2*log2(M)*10.^(EBN0dB/10))*sin(pi/M))/log2(M); % TEB QPSK Théorique
semilogy(EBN0dB,TEB,'b-')
hold on
semilogy(EBN0dB,TEB_theorique,'r-')
