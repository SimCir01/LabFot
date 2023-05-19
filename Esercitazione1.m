% ESERCITAZIONE 1
% Lo scopo di questa esercitazione è di riuscire ad ottenere, a partire
% dall'analisi degli spettri di assorbimento relativi ai dati ottenuti
% dallo spettrofotometro, una curva di calibrazione
% (Assorbanza_picco - Concentrazione) che permetta di ricavare il valore di
% concentrazione cercato per un campione di Doxorubicina.

close all;
clc;
clear all;

% Inizializzazione vettori per grafico completo delle rette di calibrazione

P = [];             % Matrice (i, j) contenente pendenze (1, j) 
                    % ed intercette (2, j) 
                    % di tutte le rette di calibrazione j-esime

Q = [];             % Vettore delle intercette delle rette di calibrazione

X_atteso = [];      % Vettore dei valori attesi di concentrazione incognita
                    % per le varie rette di calibrazione

X_lin = [];         % Vettore dei valori di concentrazione incognita
                    % ottenuti tramite calibrazione lineare per le varie
                    % rette di calibrazione

X_dil_atteso = [];  % Vettore dei valori attesi di concentrazione incognita
                    % per le varie rette di calibrazione,
                    % tenendo conto della diluizione

X_dil_lin = [];     % Vettore dei valori di concentrazione incognita
                    % ottenuti tramite calibrazione lineare per le varie
                    % rette di calibrazione, tenendo conto della diluizione

MAX = [];           % Vettore dei picchi massimi di assorbanza relativi
                    % alla concentrazione incognita per le varie
                    % rette di calibrazione

%% VIS 330-660 nm %%
%%% PLOT DEGLI SPETTRI DI ASSORBIMENTO %%%

% Import dei dati
for i=1:8
    dati = xlsread('laboratorio _050423.xlsx', i, 'C281:D941' ); 
    A(:,i) = dati(:,2);
end
lambda = dati(:,1); 


% Acquisizione degli spettri di assorbimento
figure 
for i = 1:8
    plot(lambda, A(:,i), 'LineWidth',1) 
    hold all
    xlabel('Wavelength [nm]')
    ylabel('Absorbance [A.U.]')
    title(['Spettri di assorbimento al variare delle concentrazioni:' ...
        ' VIS (330-660) nm '])
    legend('0 \mug ml', '16 \mug ml', '32 \mug ml', '40 \mug ml', ...
        '48 \mug ml','64 \mug ml', '80 \mug ml', 'campione scena crimine')
end    

% Calcolo del valore massimo di assorbanza in corrispondenza del picco 
% per ciascun valore di concentrazione
M = max(A);

%%% ELABORAZIONE DELLA RETTA DI CALIBRAZIONE %%%

% Coordinate sull'asse x dei punti sperimentali
C = [0 16 32 40 48 64 80];

% Cordinate sull'asse y dei punti sperimentali,
% ad eccezione del valore di assorbanza relativo al
% campione della scena del crimine
M_sc_VIS = M(1:7);

%Calcolo della lambda di assorbimento, alla quale si ha il picco massimo
Lambda_Amax = [];
for i = 1:8
    %indici del vettore A(:,i) dei picchi cui trovare picco max M
    index_lambda_Amax = find(A(:,i)==M(i));   
    %valori di lambda relativi a tale indice per la concentrazione i-esima
    lambda_Amax_tot = lambda(index_lambda_Amax);
    %media di questi valori
    lambda_Amax_tot = mean(lambda_Amax_tot);
    %aggiunta al vettore che memmorizza il valore medio
    Lambda_Amax(end + 1) = lambda_Amax_tot;
end
%media tra i valori medi di lamda relativi a ciascuna concentrazione
%cioè valore medio di lambda
LAMBDA_Amax = mean(Lambda_Amax);
LAMBDA_Amax = round(LAMBDA_Amax);

% Matrice 1column=lamda 2-9column=A
Mat = [];
Mat(:,1)= lambda;
for i= 2:9
    Mat(:,i) = A(i);
end

% Nota la LAMBDA_Amax, la si fissa e si calcolano i picchi ad ogni C
picchi=[];
for i=1:8
    index_lambda_Amax = find(Mat(:,1)==LAMBDA_Amax);
    picco = A(index_lambda_Amax,i);
    picco = mean(picco);
    picchi(end + 1) = picco;
end
picchi_sc_VIS = picchi(1:7);

% Calibrazione con polinomio di ordine 1 (retta)

%inizio calibrazione
figure
plot(C, picchi_sc_VIS, 'ok','LineWidth',1) 
xlabel('Concentrazione [\mug ml]')
ylabel('Absorbance [A.U.]')
title('Retta di calibrazione della Doxorubicina: VIS (330-660 nm)')
p = polyfit(C, picchi_sc_VIS, 1);
x2=[0:0.1:80];
y2=polyval(p, x2);
hold on
plot(x2, y2, 'k', 'LineWidth',1) 
%fine calibrazione

%%% CALCOLO DELLA CONCENTRAZIONE DEL CAMPIONE %%%

% valore di concentrazione calcolato seguendo la formula teorica per
% l'interpolazione
x_atteso = ((picchi(8)-picchi(5))*(C(6)-C(5)))/(picchi(6)-picchi(5)) + C(5);


% pendenza della retta di calibrazione 
m = p(1);
% intercetta all'origine della retta di calibrazione
q = p(2);

% valore di concentrazione del campione della scena del crimine
x_lin = (picchi(8)-q)/m;
plot(x_lin, picchi(8), 'or','LineWidth', 1);
plot(x_atteso, picchi(8), '*r','LineWidth',1);
legend('Punti sperimentali', 'Linear', ...
    'Incognita con calibrazione lineare', 'Incognita Attesa')
grid on

% valore di concentrazione del campione incognito tenuto conto della
% diluizione (5000 volte)

x_dil_atteso = 5000*x_atteso;
x_dil_lin = 5000*x_lin;

P(:, end + 1) = p;
Q(:, end + 1) = q;
X_atteso(:, end + 1) = x_atteso;
X_lin(:, end + 1) = x_lin;
MAX(:, end + 1) = M(8);

%% UV-VIS 190-850 nm %%

%%% PLOT DEGLI SPETTRI DI ASSORBIMENTO %%%

A = [];
% Import dei dati
for i=1:8
    dati = xlsread('laboratorio _050423.xlsx', i, 'C1:D1321' ); 
    A(:,i) = dati(:,2);
end
lambda = dati(:,1); 

% Acquisizione degli spettri di assorbimento
figure 
for i = 1:8
    plot(lambda, A(:,i),'LineWidth',1)  
    hold all
    xlabel('Wavelength [nm]')
    ylabel('Absorbance [A.U.]')
    title(['Spettri di assorbimento al variare delle concentrazioni:' ...
        ' UV-VIS (190-850 nm)'])
    legend('0 \mug ml', '16 \mug ml', '32 \mug ml', '40 \mug ml', ...
        '48 \mug ml','64 \mug ml', '80 \mug ml', 'campione scena crimine')
end    

% Calcolo del valore massimo di assorbanza in corrispondenza del picco 
% per ciascun valore di concentrazione
M = max(A);


%%% ELABORAZIONE DELLA RETTA DI CALIBRAZIONE %%%

% Coordinate sull'asse x dei punti sperimentali
C = [0 16 32 40 48 64 80];

% Cordinate sull'asse y dei punti sperimentali,
% ad eccezione del valore di assorbanza relativo al
% campione della scena del crimine
M_sc_UV_VIS = M(1:7);


% Calibrazione con polinomio di ordine 1 (retta)

%inizio calibrazione
figure
plot(C, M_sc_UV_VIS, 'ok','LineWidth',1) 
xlabel('Concentrazione [\mug ml]')
ylabel('Absorbance [A.U.]')
title('Retta di calibrazione della Doxorubicina: UV-VIS (190-850 nm)')
p = polyfit(C, M_sc_UV_VIS, 1);
x2=[0:0.1:80];
y2=polyval(p, x2);
hold on
plot(x2, y2, 'k','LineWidth',1) 
%fine calibrazione

%%% CALCOLO DELLA CONCENTRAZIONE DEL CAMPIONE %%%

% valore di concentrazione calcolato seguendo la formula teorica per
% l'interpolazione
x_atteso = ((M(8)-M(5))*(C(6)-C(5)))/(M(6)-M(5)) + C(5);


% pendenza della retta di calibrazione 
m = p(1);
% intercetta all'origine della retta di calibrazione
q = p(2);

% valore di concentrazione del campione della scena del crimine
x_lin = (M(8)-q)/m;
plot(x_lin, M(8), 'or','LineWidth',1) ;
plot(x_atteso, M(8), '*r','LineWidth',1) 
legend('Punti sperimentali', 'Linear', ...
    'Incognita con calibrazione lineare', 'Incognita Attesa')
grid on

% valore di concentrazione del campione incognito tenuto conto della
% diluizione (5000 volte)

x_dil_atteso = 5000*x_atteso;
x_dil_lin = 5000*x_lin;

P(:, end + 1) = p;
Q(:, end + 1) = q;
X_atteso(:, end + 1) = x_atteso;
X_lin(:, end + 1) = x_lin;
MAX(:, end + 1) = M(8);

%% UV 284-325 nm %%

%%% PLOT DEGLI SPETTRI DI ASSORBIMENTO %%%
A = [];
% Import dei dati
for i=1:8
    dati = xlsread('laboratorio _050423.xlsx', i, 'C189:D271' ); 
    A(:,i) = dati(:,2);
end
lambda = dati(:,1); 

% Acquisizione degli spettri di assorbimento
figure 
for i = 1:8
    plot(lambda, A(:,i), 'LineWidth',1)  
    hold all
    xlabel('Wavelength [nm]')
    ylabel('Absorbance [A.U.]')
    title(['Spettri di assorbimento al variare delle concentrazioni:' ...
        ' UV (284-325) nm'])
    legend('0 \mug ml', '16 \mug ml', '32 \mug ml', '40 \mug ml', ...
        '48 \mug ml','64 \mug ml', '80 \mug ml', 'campione scena crimine')
end    

% Calcolo del valore massimo di assorbanza in corrispondenza del picco 
% per ciascun valore di concentrazione
M = max(A);


%%% ELABORAZIONE DELLA RETTA DI CALIBRAZIONE %%%

% Coordinate sull'asse x dei punti sperimentali
C = [0 16 32 40 48 64 80];

% Cordinate sull'asse y dei punti sperimentali,
% ad eccezione del valore di assorbanza relativo al
% campione della scena del crimine
M_sc_UV = M(1:7);

%Calcolo della lambda di assorbimento, alla quale si ha il picco massimo
Lambda_Amax = [];
for i = 1:8
    %indici del vettore A(:,i) dei picchi cui trovare picco max M
    index_lambda_Amax = find(A(:,i)==M(i));   
    %valori di lambda relativi a tale indice per la concentrazione i-esima
    lambda_Amax_tot = lambda(index_lambda_Amax);
    %media di questi valori
    lambda_Amax_tot = mean(lambda_Amax_tot);
    %aggiunta al vettore che memmorizza il valore medio
    Lambda_Amax(end + 1) = lambda_Amax_tot;
end
%media tra i valori medi di lamda relativi a ciascuna concentrazione
%cioè valore medio di lambda
LAMBDA_Amax = mean(Lambda_Amax);
LAMBDA_Amax = round(LAMBDA_Amax);

% Matrice 1column=lamda 2-9column=A
Mat = [];
Mat(:,1)= lambda;
for i= 2:9
    Mat(:,i) = A(i);
end

% Nota la LAMBDA_Amax, la si fissa e si calcolano i picchi ad ogni C
picchi=[];
for i=1:8
    index_lambda_Amax = find(Mat(:,1)==LAMBDA_Amax);
    picco = A(index_lambda_Amax,i);
    picco = mean(picco);
    picchi(end + 1) = picco;
end
picchi_sc_UV = picchi(1:7);


% Calibrazione con polinomio di ordine 1 (retta)

%inizio calibrazione
figure
plot(C, picchi_sc_UV, 'ok', 'LineWidth',1) 
xlabel('Concentrazione [\mug ml]')
ylabel('Absorbance [A.U.]')
title(['Retta di calibrazione della Doxorubicina: ' ...
    'UV (284-325 nm)'])
p = polyfit(C, picchi_sc_UV, 1);
x2=[0:0.1:80];
y2=polyval(p, x2);
hold on
plot(x2, y2, 'k', 'LineWidth',1) 
%fine calibrazione

%%% CALCOLO DELLA CONCENTRAZIONE DEL CAMPIONE %%%

% valore di concentrazione calcolato seguendo la formula teorica per
% l'interpolazione
x_atteso = ((picchi(8)-picchi(5))*(C(6)-C(5)))/(picchi(6)-picchi(5))+C(5);


% pendenza della retta di calibrazione 
m = p(1);
% intercetta all'origine della retta di calibrazione
q = p(2);

% valore di concentrazione del campione della scena del crimine
x_lin = (picchi(8)-q)/m;
plot(x_lin, picchi(8), 'or', 'LineWidth',1) ;
plot(x_atteso, picchi(8), '*r', 'LineWidth',1);
legend('Punti sperimentali', 'Linear', ...
    'Incognita con calibrazione lineare', 'Incognita Attesa')
grid on

% valore di concentrazione del campione incognito tenuto conto della
% diluizione (5000 volte)

x_dil_atteso = 5000*x_atteso;
x_dil_lin = 5000*x_lin;

P(:, end + 1) = p;
Q(:, end + 1) = q;
X_atteso(:, end + 1) = x_atteso;
X_lin(:, end + 1) = x_lin;
MAX(:, end + 1) = M(8);

%% PLOT DELLE RETTE DI CALIBRAZIONE DEI DIVERSI RANGE DI WAVELENGTH

figure 
x = 0:0.1:80;
Y1 = P(1,1)*x+P(2,1);                   % retta di calibrazione VIS
Y2 = P(1,2)*x+P(2,2);                   % retta di calibrazione UV-VIS
Y3 = P(1,3)*x+P(2,3);                   % retta di calibrazione UV



% Plot della retta di calibrazione VIS

p1_lin = plot(X_lin(1), MAX(1), 'xr', 'LineWidth',1);
hold on
plot(C, picchi_sc_VIS, 'xb','LineWidth',1)
%plot(X_atteso(1), MAX(1), '*b');
%hold on
p1 = plot(x, Y1, '-b', 'LineWidth',1);
p_VIS = polyfit(C, picchi_sc_UV, 1);
x2=[0:0.1:80];
y2_VIS=polyval(p_VIS, x2);
hold on
plot(x2, y2_VIS, '-b','LineWidth',1)
hold on


% Plot della retta di calibrazione UV-VIS

% p2_lin = plot(X_lin(2), MAX(2), 'om', 'LineWidth',1) ;
% hold on
% plot(C, M_sc_UV_VIS, 'om','LineWidth',1)
% %plot(X_atteso(2), MAX(2), '*m');
% %hold on
% p2 = plot(x, Y2, '-m', 'LineWidth',1) ;
% p_UV_VIS = polyfit(C, M_sc_UV, 1);
% x2=[0:0.1:80];
% y2_UV_VIS=polyval(p_UV_VIS, x2);
% hold on
% plot(x2, y2_UV_VIS, '-m','LineWidth',1)
% hold on


% Plot della retta di calibrazione UV

p3_lin = plot(X_lin(3), MAX(3), 'or', 'LineWidth',1) ;
hold on
plot(C, picchi_sc_UV, 'og','LineWidth',1)
%plot(X_atteso(3), MAX(3), '*g');
%hold on
p3 = plot(x, Y3, '-g', 'LineWidth',1);
p_UV = polyfit(C, picchi_sc_UV, 1);
x2=[0:0.1:80];
y2_UV=polyval(p_UV, x2);
hold on
plot(x2, y2_UV, '-g','LineWidth',1)


%h = [p1;p2;p3;p1_lin;p2_lin;p3_lin];
h = [p1;p3;p1_lin;p3_lin];
s_X_lin1 = sprintf('Concentrazione incognita VIS: %.2f', X_lin(1));
%s_X_lin2 = sprintf('Concentrazione incognita: %.2f', X_lin(2));
s_X_lin3 = sprintf('Concentrazione incognita UV: %.2f', X_lin(3));
%legend(h, 'VIS','UV-VIS', 'UV', s_X_lin1, s_X_lin2, s_X_lin3)
legend(h, 'VIS', 'UV', s_X_lin1, s_X_lin3)
title('Rette di calibrazione')

