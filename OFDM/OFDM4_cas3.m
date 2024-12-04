close all;
clear all;

% Variables d'initialisation

N = 16;     % Nombre de porteuses   
Nbits = 1000;   % Nombre de bits total
Nactives = 16;  % Nombre de porteuses actives
lignes_garde = 6;   % Nombre de lignes de garde

% Module & phase de la réponse en fréquence du canal de propagation

H = [0.407,0.815,0.407];
figure('Name','Module & phase de la réponse en fréquence du canal de propagation')
freqz(H,1,1024,"whole")

% Mapping du signal

Mapping = zeros(N,Nbits);
for i=1:Nactives
    Mapping(i,:) = randi([0 1],1,Nbits)*2 - 1;  % BPSK
end

% Canal Proakis

Xe = ifft(Mapping,N);

garde = Xe(N-lignes_garde+1:end,:); 
Xe_garde = [garde;Xe];  % On introduit ici un préfixe cyclique "au dessus" de Xe

Xe_temp = reshape(Xe_garde,1,(N+lignes_garde)*Nbits);   % On met la matrice Xe + l'intervalle de garde en ligne
Xe_temp(1:7) = [];                                      % On supprime 7 éléments au début
Xe_temp((N+lignes_garde)*Nbits-21:(N+lignes_garde)*Nbits-7)=[];     % On en supprime 15 à la fin : au total un symbole a été supprimé (22 éléments)

Xe_cas3 = reshape(Xe_temp,22,Nbits-1);              % On reshape la matrice avec le dernier symbole supprimé

ofdm_lineaire = reshape(Xe_cas3, 1, (Nbits-1)*(N+lignes_garde));
SignalSortieCanal=filter(H,1,ofdm_lineaire);


% Figure DSP

DSP = pwelch(SignalSortieCanal,[],[],[],16,'centered');

figure('Name','DSP')
plot(10*log(DSP))
xlabel('Fréquence (Hz)')
ylabel('DSP')

% Calcul de Ck

C_k = fft(H,16);                % Calcul des coefficients C(k)
Egalisateur_ML = repmat(C_k(:),1,Nbits-1);   % Ajustement de la taille de la matrice

% Démodulation

Y_reshape = reshape(SignalSortieCanal,size(Xe_cas3));      % On reshape la matrice en ligne
Y_recu = fft(Y_reshape,N);
Y_egalisation = conj(Egalisateur_ML).*Y_recu;       % On multiplie les éléments de Y_recu par les éléments correspondants de 1/H_k
                                        % Pour contrer les effets du canal

% Constellations 6 et 15

FPorteuse6 = Y_egalisation(6,:);   % On extrait les points sur la porteuse 6
FPorteuse15 = Y_egalisation(15,:); % On extrait les points sur la porteuse 15

figure('Name','Constellation de la 6ème porteuse');
scatter(real(FPorteuse6), imag(FPorteuse6));    % On affiche la partie imaginaire en fonction de la partie réelle de chaque point
xlabel('Partie réelle');
ylabel('Partie imaginaire');

figure('Name','Constellation de la 15ème porteuse');
scatter(real(FPorteuse15), imag(FPorteuse15));  % On affiche la partie imaginaire en fonction de la partie réelle de chaque point
xlabel('Partie réelle');
ylabel('Partie imaginaire');

