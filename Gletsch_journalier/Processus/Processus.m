%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%/Ce programme simule l'effet des trois fonction de transfert obtenu par décomposition
% en élément simple de la fonction suivante H(z)=    0.051
%                                                 ---------------- 
%                                                (z+0.9655)(z-0.05299)
%On obtient : H(z)= K + H1 + H2  (Interpréter comme 3 processus)
%//////////////////////////////////////////////////////////////////////////
%/// H(Z) est la fonction de transfert d'un modèle ARX(2,1) obtenu par/////
%/// la méthode MMC sur la crue de l'année 1999 ///////////////////////////
%//////////////////////////////////////////////////////////////////////////
%**************************************************************************


function Processus

t=[1:100];      %// Echelle de temps : 100 pas
e=(exp(-(t-10).*(t-10))/10);  %// Entrée, Exitation 


%sys=tf([0.05],[1,0.9655],1)    %// Fonction de transfert H1 de l'écoulement lent
sys=tf([7.4],[1,21])

%sys2=tf([0.9578],[1,-0.0529],1)  %// Fonction de Transfert H2 de l'écoulement rapide
sys2=tf([-0.32])

%sys3=tf([-0.028],[1],1)         %// Interpréter K comme l'effet d'absorption, infiltration

figure;
lsim(sys,e,t)  %// Simulation du processus d'écoulement lent
hold on,lsim(sys2,e,t)  %// Simulation du procesus d'écouolement rapide
%hold on,lsim(sys3,e)  %// Simulation de l'effet d'infiltration