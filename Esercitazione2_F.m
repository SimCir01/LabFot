close all;
clc;
clear all;


%% MISURE DI FLUORESCENZA
%%% PLOT DEGLI SPETTRI DI ASSORBIMENTO %%%
lambda = xlsread('FLORESCENZE.xlsx', 'A7:A1618' );
dark = xlsread('FLORESCENZE.xlsx', 'B7:B1618' );
C5_F = xlsread('FLORESCENZE.xlsx', 'C7:C1618' );
C10_F = xlsread('FLORESCENZE.xlsx', 'D7:D1618' );
C20_F = xlsread('FLORESCENZE.xlsx', 'E7:E1618' );
C30_F = xlsread('FLORESCENZE.xlsx', 'F7:F1618' );
CX_F = xlsread('FLORESCENZE.xlsx', 'G7:G1618' );
CX1_F = xlsread('FLORESCENZE.xlsx', 'H7:H1618' );
CX2_F = xlsread('FLORESCENZE.xlsx', 'I7:I1618' );

figure 
hold all
plot(lambda,C5_F)
plot(lambda, C10_F)
plot(lambda, C20_F)
plot(lambda, C30_F)
plot(lambda, CX_F)
plot(lambda, CX1_F)
plot(lambda, CX2_F)
%axis([175 1200 -5 10])
xlabel('Wavelength [nm]')
ylabel('Absorbance [A.U.]')
title('Spettri di fluorescenza al variare delle concentrazioni (dark non sottratto)')
legend('5 \mug ml', '1O \mug ml', '20 \mug ml', '30 \mug ml', ...
        'X \mug ml','X1 \mug ml', 'X2 \mug ml')

% Sottrazione del dark
C5_F = C5_F - dark;
C10_F = C10_F - dark;
C20_F = C20_F - dark;
C30_F = C30_F - dark;
CX_F = CX_F - dark;
CX1_F = CX1_F - dark;
CX2_F = CX2_F - dark;

% Media tra le concentrazioni incognite
CXM_F = (CX_F + CX1_F + CX2_F)/3;

%Acquisizione degli spettri di assorbimento
figure 
hold all
plot(lambda,C5_F)
plot(lambda, C10_F)
plot(lambda, C20_F)
plot(lambda, C30_F)
plot(lambda, CX_F)
plot(lambda, CX1_F)
plot(lambda, CX2_F)
%axis([175 1200 -5 10])
xlabel('Wavelength [nm]')
ylabel('Absorbance [A.U.]')
title('Spettri di fluorescenza al variare delle concentrazioni (dark sottratto) ')
legend('5 \mug ml', '10 \mug ml', '20 \mug ml', '30 \mug ml', ...
        'X \mug ml','X1 \mug ml', 'X2 \mug ml')
  
P = [];
X_atteso = [];
X_lin = [];
CI=[];
sign=['x','o', '^'];
j = 1;
C_inc = [CX_F, CX1_F, CX2_F];
figure
hold on
for i=1:length(C_inc(1,:))
    %% CX_A
    % Matrice A dei picchi di assorbanza per ciascuna C
    A=[C5_F, C10_F, C20_F, C30_F, C_inc(:,i)];
    
    % Calcolo dei massimi
    M = max(A);
    
    %%% ELABORAZIONE DELLA RETTA DI CALIBRAZIONE %%%
    
    % Coordinate sull'asse x dei punti sperimentali
    C = [5 10 20 30];
    
    % Cordinate sull'asse y dei punti sperimentali,
    % ad eccezione del valore di assorbanza relativo al
    % campione incognito
    M_sc = M(1:4);
    
%     %Calcolo della lambda di assorbimento, alla quale si ha il picco massimo
%     Lambda_Amax = [];
%     for i=1:5
%     %indici del vettore A(:,i) dei picchi cui trovare picco max M
%     index_lambda_Amax = find(A(:,i)==M(i));
%     %valori di lambda relativi a tale indice per la concentrazione i-esima
%     lambda_Amax_tot = lambda(index_lambda_Amax);
%     %aggiunta al vettore che memmorizza il valore medio
%     Lambda_Amax(end + 1) = lambda_Amax_tot;
%     end
% 
% %     
%     % Media tra i valori medi di lamda relativi a ciascuna concentrazione
%     %cioè valore medio di lambda
%     LAMBDA_Amax = mean(Lambda_Amax);
%     
% %     
%     % Matrice 1column=lambda 2-7column=A
%     Mat = [];
%     Mat(:,1)= lambda;
%     for i= 2:6
%         Mat(:,i) = A(i-1);
%     end
% %     
%     % Nota la LAMBDA_Amax, la si fissa e si calcolano i picchi ad ogni C
%     picchi=[];
%     for i=1:1
%         index_lambda_Amax = find(Mat(:,1)==LAMBDA_Amax);
%         picco = A(index_lambda_Amax,i);
%         picco = mean(picco);
%         picchi(end + 1) = picco;
%     end
%     picchi_sc = picchi(1:4);
    
    % Calibrazione con polinomio di ordine 1 (retta)
    
    %inizio calibrazione
    
    plot(C, M_sc, 'ok','LineWidth',1) 
    xlabel('Concentrazione [\mug ml]')
    ylabel('Fluorescence [A.U.]')
    title('Retta di calibrazione della Rodamina')
    p = polyfit(C, M_sc, 1);
    x2=[0:0.1:35];
    y2=polyval(p, x2);
    hold on
    plot(x2, y2, '-k','LineWidth',1 );
    hold on
    %fine calibrazione
    
    %%% CALCOLO DELLA CONCENTRAZIONE DEL CAMPIONE %%%
    
    % valore di concentrazione calcolato seguendo la formula teorica per
    % l'interpolazione
    x_atteso = ((M(5)-M(2))*(C(3)-C(2)))/(M(3)-M(2)) + C(2);
    X_atteso(end+1) = x_atteso;
    
    % pendenza della retta di calibrazione 
    m = p(1);
    % intercetta all'origine della retta di calibrazione
    q = p(2);
    
    % valore di concentrazione del campione della scena del crimine
    x_lin = (M(5)-q)/m;
    X_lin(end+1) = x_lin;
    CI(end+1) = plot(x_lin, M(5), [sign(j), 'r']);
    %plot(x_atteso, picchi(5), '*r','LineWidth',1);
    grid on
    j=j+1;
end
s_X_lin1 = sprintf('Concentrazione incognita AX: %.2f', X_lin(1));
s_X_lin2 = sprintf('Concentrazione incognita AX1: %.2f', X_lin(2));
s_X_lin3 = sprintf('Concentrazione incognita AX2: %.2f', X_lin(3));
h = [CI(1), CI(2), CI(3)];
legend(h, s_X_lin1, s_X_lin2, s_X_lin3)

% Calcolo della media e dell'incertezza della concentrazione incognita
somma = 0;
X_lin_m = mean(X_lin);
for z = 1:3
    iter = (X_lin(z) - X_lin_m)^2;
    somma = somma + iter;
end

scarto_t = sqrt(0.5*somma);
u = scarto_t/(sqrt(3));

% Si assume distribuzione uniforme (rettangolare) con L.C. del 90%
LC = 90;
k = (LC/100)*sqrt(3);
U = k*u;

