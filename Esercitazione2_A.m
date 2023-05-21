% ESERCITAZIONE 2

close all;
clc;
clear all;

%% MISURE DI ASSORBANZA
%%% PLOT DEGLI SPETTRI DI ASSORBIMENTO %%%
lambda = xlsread('ASSORBANZE.xlsx', 'A7:A1618' ); 
C5_A = xlsread('ASSORBANZE.xlsx', 'B7:B1618' );
C10_A = xlsread('ASSORBANZE.xlsx', 'C7:C1618' );
C20_A = xlsread('ASSORBANZE.xlsx', 'D7:D1618' );
C30_A = xlsread('ASSORBANZE.xlsx', 'E7:E1618' );
CX_A = xlsread('ASSORBANZE.xlsx', 'F7:F1618' );
CX1_A = xlsread('ASSORBANZE.xlsx', 'G7:G1618' );
CX2_A = xlsread('ASSORBANZE.xlsx', 'H7:H1618' );
CXM_A = (CX_A + CX1_A + CX2_A)/3;

%Acquisizione degli spettri di assorbimento
figure 
hold all
plot(lambda,C5_A)
plot(lambda, C10_A)
plot(lambda, C20_A)
plot(lambda, C30_A)
plot(lambda, CX_A)
plot(lambda, CX1_A)
plot(lambda, CX2_A)
%axis([175 1200 -5 10])
xlabel('Wavelength [nm]')
ylabel('Absorbance [A.U.]')
title('Spettri di assorbimento al variare delle concentrazioni')
legend('5 \mug ml', '1O \mug ml', '20 \mug ml', '30 \mug ml', ...
        'X \mug ml','X1 \mug ml', 'X2 \mug ml')

% Range di lunghezze d'onda da considerare per l'assorbanza
lambda = lambda(237:338);

% Picchi di A rispettivi al range di lambda considerato
C5_A = C5_A(237:338);
C10_A = C10_A(237:338);
C20_A = C20_A(237:338);
C30_A = C30_A(237:338);
CX_A = CX_A(237:338);
CX1_A = CX1_A(237:338);
CX2_A = CX2_A(237:338);
CXM_A = CXM_A(237:338);

% figure 
% hold all
% plot(lambda,C5_A)
% plot(lambda, C10_A)
% plot(lambda, C20_A)
% plot(lambda, C30_A)
% plot(lambda, CX_A)
% plot(lambda, CX1_A)
% plot(lambda, CX2_A)
% %axis([175 1200 -5 10])
% xlabel('Wavelength [nm]')
% ylabel('Absorbance [A.U.]')
% title('Spettri di assorbimento al variare delle concentrazioni')
% legend('5 \mug ml', '1O \mug ml', '20 \mug ml', '30 \mug ml', ...
%         'X \mug ml','X1 \mug ml', 'X2 \mug ml')

P = [];
X_atteso = [];
X_lin = [];
CI=[];
sign=['x','o', '^'];
j = 1;
C_inc = [CX_A, CX1_A, CX2_A];
figure
hold on
for i=1:length(C_inc(1,:))
    %% CX_A
    % Matrice A dei picchi di assorbanza per ciascuna C
    A=[C5_A, C10_A, C20_A, C30_A, C_inc(:,i)];
    
    % Calcolo dei massimi
    M = max(A);
    
    %%% ELABORAZIONE DELLA RETTA DI CALIBRAZIONE %%%
    
    % Coordinate sull'asse x dei punti sperimentali
    C = [5 10 20 30];
    
    % Cordinate sull'asse y dei punti sperimentali,
    % ad eccezione del valore di assorbanza relativo al
    % campione incognito
    M_sc = M(1:4);
    
    %Calcolo della lambda di assorbimento, alla quale si ha il picco massimo
    Lambda_Amax = [];
    for i=1:5
    %indici del vettore A(:,i) dei picchi cui trovare picco max M
    index_lambda_Amax = find(A(:,i)==M(i));
    %valori di lambda relativi a tale indice per la concentrazione i-esima
    lambda_Amax_tot = lambda(index_lambda_Amax);
    %aggiunta al vettore che memmorizza il valore medio
    Lambda_Amax(end + 1) = lambda_Amax_tot;
    end
    
    % Media tra i valori medi di lamda relativi a ciascuna concentrazione
    %cioè valore medio di lambda
    LAMBDA_Amax = mean(Lambda_Amax);
    
    % Matrice 1column=lambda 2-7column=A
    Mat = [];
    Mat(:,1)= lambda;
    for i= 2:6
        Mat(:,i) = A(i-1);
    end
    
    % Nota la LAMBDA_Amax, la si fissa e si calcolano i picchi ad ogni C
    picchi=[];
    for i=1:5
        index_lambda_Amax = find(Mat(:,1)==LAMBDA_Amax);
        picco = A(index_lambda_Amax,i);
        %picco = mean(picco);
        picchi(end + 1) = picco;
    end
    picchi_sc = picchi(1:4);
    
    % Calibrazione con polinomio di ordine 1 (retta)
    
    %inizio calibrazione
    
    plot(C, picchi_sc, 'ok','LineWidth',1) 
    xlabel('Concentrazione [\mug ml]')
    ylabel('Absorbance [A.U.]')
    title('Retta di calibrazione della Rodamina')
    p = polyfit(C, picchi_sc, 1);
    x2=[0:0.1:35];
    y2=polyval(p, x2);
    hold on
    plot(x2, y2, '-k','LineWidth',1 );
    hold on
    %fine calibrazione
    
    %%% CALCOLO DELLA CONCENTRAZIONE DEL CAMPIONE %%%
    
    % valore di concentrazione calcolato seguendo la formula teorica per
    % l'interpolazione
    x_atteso = ((picchi(5)-picchi(2))*(C(3)-C(2)))/(picchi(3)-picchi(2)) + C(2);
    X_atteso(end+1) = x_atteso;
    
    % pendenza della retta di calibrazione 
    m = p(1);
    % intercetta all'origine della retta di calibrazione
    q = p(2);
    
    % valore di concentrazione del campione della scena del crimine
    x_lin = (picchi(5)-q)/m;
    X_lin(end+1) = x_lin;
    CI(end+1) = plot(x_lin, picchi(5), [sign(j), 'r']);
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



