%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%/Ce programme simule l'effet des trois fonction de transfert obtenu par d�composition
% en �l�ment simple de la fonction suivante H(z)=    0.051
%                                                 ---------------- 
%                                                (z+0.9655)(z-0.05299)
%On obtient : H(z)= K + H1 + H2  (Interpr�ter comme 3 processus)
%//////////////////////////////////////////////////////////////////////////
%/// H(Z) est la fonction de transfert d'un mod�le ARX(2,1) obtenu par/////
%/// la m�thode MMC sur la crue de l'ann�e 1999 ///////////////////////////
%//////////////////////////////////////////////////////////////////////////
%**************************************************************************


function Processus

t=[1:100];      %// Echelle de temps : 100 pas
e=(exp(-(t-10).*(t-10))/10);  %// Entr�e, Exitation 


%sys=tf([0.05],[1,0.9655],1)    %// Fonction de transfert H1 de l'�coulement lent
sys=tf([7.4],[1,21])

%sys2=tf([0.9578],[1,-0.0529],1)  %// Fonction de Transfert H2 de l'�coulement rapide
sys2=tf([-0.32])

%sys3=tf([-0.028],[1],1)         %// Interpr�ter K comme l'effet d'absorption, infiltration

figure;
lsim(sys,e,t)  %// Simulation du processus d'�coulement lent
hold on,lsim(sys2,e,t)  %// Simulation du procesus d'�couolement rapide
%hold on,lsim(sys3,e)  %// Simulation de l'effet d'infiltration