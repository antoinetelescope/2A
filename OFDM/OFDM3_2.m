close all;
clear all;

% Variables d'initialisation

N = 16;     % Nombre de porteuses   
Nbits = 1000;   % Nombre de bits total
Nactives = 16;  % Nombre de porteuses actives
lignes_garde = 2;   % Nombre de lignes de garde

% Module & phase de la réponse en fréquence du canal de propagation

H = [0.407,0.815,0.407];
figure('Name','Module & phase de la réponse en fréquence du canal de propagation')
freqz(H,1,1024,"whole")

% Mapping du signal

Mapping = zeros(N,Nbits);
for i=1:Nactives
    Mapping(i,:) = randi([0 1],1,Nbits)*2 - 1;  %BPSK
end

% Canal Proakis

Xe = ifft(Mapping,N);

garde = zeros(lignes_garde, Nbits); 
Xe_garde = [garde;Xe];  % On introduit ici 2 lignes de zéro au dessus de Xe : intervalle de garde


ofdm_lineaire = reshape(Xe_garde, 1, Nbits*(N+lignes_garde));
SignalSortieCanal=filter(H,1,ofdm_lineaire);


% Figure DSP

DSP = pwelch(SignalSortieCanal,[],[],[],16,'centered');

figure('Name','DSP')
plot(10*log(DSP))
xlabel('Fréquence (Hz)')
ylabel('DSP')

% Démodulation

Y_reshape = reshape(SignalSortieCanal,size(Xe_garde));      % On reshape la matrice en ligne
Xe_sansgarde = Y_reshape(lignes_garde+1:N+lignes_garde,:);      % On supprime les 2 lignes de zéro
Y_recu = fft(Xe_sansgarde,N);

% Constellations 6 et 15

FPorteuse6 = Y_recu(6,:);   % On extrait les points sur la porteuse 6
FPorteuse15 = Y_recu(15,:); % On extrait les points sur la porteuse 15

figure('Name','Constellation de la 6ème porteuse');
scatter(real(FPorteuse6), imag(FPorteuse6));    % On affiche la partie imaginaire en fonction de la partie réelle de chaque point
xlabel('Partie réelle');
ylabel('Partie imaginaire');

figure('Name','Constellation de la 15ème porteuse');
scatter(real(FPorteuse15), imag(FPorteuse15));  % On affiche la partie imaginaire en fonction de la partie réelle de chaque point
xlabel('Partie réelle');
ylabel('Partie imaginaire');

% TEB

Y_fin = real(Y_recu)>0;
Y_fin = Y_fin*2-1;
TEB = mean(Mapping~=Y_fin,"all");

disp('Le TEB est de');
disp(TEB);


