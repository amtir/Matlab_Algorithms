



%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%///// A.A.P : Algorithme d'adaptation paramétrique ///////////////////////
%//////////////////////////////////////////////////////////////////////////
%///// Methode des moindres carres ponderes recurrents/////////////////////
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%///// FILTRE NUMERIQUE NON-STATIONAIRE, paramètres adaptatifs ////////////
%///// FILTRE A MEMOIRE LIMITEE ///////////////////////////////////////////
%///// Estimation des paramètres du modèle ARX :  A(z)Y(k)=B(Z)U(k)+e(k) //
%///// une variable explicative U=P(précipitation)  ///////////////////////
%///// la sortie est le débit Q, d'ordre ARX(na,nb,nk) ////////////////////
%///// A(z)Y(k)=B(Z)U(k)+e(k) /////////////////////////////////////////////
%///// Periode de prevision N1=[N1(1),N1(2)]///////////////////////////////
%///// na et nb l'ordre du model //////////////////////////////////////////
%///// nk est le retard de la variable explictive  ////////////////////////
%///// 0.95<Lamda<0.99 est le facteur d'oubli /////////////////////////////
%///// Exemple na=3, nb=3, nk=0 et N1=[2000,3500]  ////////////////////////
%//////////////////////////////////////////////////////////////////////////
% REFERENCE: R Longchamp, Commande Numérique de Systèmes Dynamiques
%                         Presse Polytechniques et universitaires Romandes
%**************************************************************************


function MMCRP


%**************************************************************************
%////////// Programme principale ///////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%////////  1 Sous fonction : Charger  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[P,Q]=Charger;

%//////////////////////////////////////////////////////////////////////////
%////////  2 Sous fonction : Partie interactive  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
sprintf('Ce programme permet d estimer le vecteur paramètre (A.A.P: MMCRP) et de simuler le débit') 
reply = input('Voulez-Vous continuer? oui=1, non=0 :');
if isempty(reply)
    quit;
elseif reply
N1 = input('Choisir la période de crue pour la simulation? : N1=[N1(1),N1(2)]=');
na = input('Choisir l ordre pour AR, na=');
nb = input('Choisir l ordre pour la partie explicative, nb=');
nk = input('Choisir le retard, si vous ne savez pas, mettre nk=0, nk=');
Lamda = input('Choisir le coefficient de pondération entre [0.98,0.99], Lamda=');
H = waitbar(0,'SVP, Attendre la fin des calculs...');


%[Q]=Methode_Lissage(Q,20,1,0.9);   %// Filtrage du debit 

if (na==0&nb==0)
    sprintf('Erreur!!! Choisir l Ordre: na, nb >0.')   
elseif (na>0 & nb==0)  %// Modele AR: Ay(t)=e(t)
     sprintf('Moddele AR') 
elseif (na==0 & nb>0)  %// y(t)=Bu(t)+e(t) 
      sprintf('Moddele X')    
elseif na>0|nb>0  %// Modele ARMAX,ARX, AR, MAX, X
        sprintf('Moddele ARX') 
end
%//////////////////////////////////////////////////////////////////////////
%////////  3 Sous fonction : Amorce de l'algorithme  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[Pm,Teta]=ConditionInitiale(na,nb,nk,N1,Q,P);

%//////////////////////////////////////////////////////////////////////////
%////////  4 Sous fonction : Algorithme de mise à jour  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[yh,A]=AlgorMiseAJour(N1,P,Q,na,nb,nk,Pm,Teta,Lamda);

%//////////////////////////////////////////////////////////////////////////
%////////  5 Sous fonction : Dessiner  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
Dessiner(yh,Q,N1,na,nb,A);


%//////////////////////////////////////////////////////////////////////////
%////////  6 Sous fonction : Criter_Validation  ///////////////////////////
%//////////////////////////////////////////////////////////////////////////
Critere_Validation(Q,yh,N1,na,nb)

%//////////////////////////////////////////////////////////////////////////
%////////  7 Sous fonction : Modele  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%Modele(N1,A,na,nb,Q,P,nk);

close(H)
end
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*************************************************************************
%//////////////////////////////////////////////////////////////////////////
%Commencer par charger les données: moyenne journalière du bassin /////////
function [P,Q]=Charger


load('Peq_nonglacier_gletsch_horaire_90_95.txt');  
Peq_nonglacier_gletsch=Peq_nonglacier_gletsch_horaire_90_95; % (32000:end);

%// Charger les Pluie equivalente de la partie glaciere (pluie +fonte neige)
load('Peq_glacier_gletsch_horaire_90_95.txt');  
Peq_glacier_gletsch=Peq_glacier_gletsch_horaire_90_95; % (32000:end);

P=(Peq_nonglacier_gletsch+Peq_glacier_gletsch);

load('Qobs_gletsch_horaire_90_95.txt'); %Débits à l'exutoire des deux bassins (Gletsh/Alptal)
Q=Qobs_gletsch_horaire_90_95; %(32000:end);  %Débits du Bassin (1)Alptal, (2)Gletsch

figure,subplot(2,1,1), plot(Q), 
ylabel('Debits à l exutoire du Bassin Gletsch (mm)'), grid on
subplot(2,1,2), plot(P)
ylabel('Precipitation Equivalente (mm)'), grid on
xlabel('Temps en pas d1 heure')%//////////////

%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************

%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
%//// Amorçage de l'algorithme à partir de N2  ////////////////////////////
%//////////////////////////////////////////////////////////////////////////

function [Pm,Teta]=ConditionInitiale(na,nb,nk,N1,Q,P);


%L=max(na,(nb+nk))+1;  %// Utiliser le maximum d'historique
L=N1(1)-(400*24);          %// Prendre une année (400jour) d'historique
%L=100;                 %// Commencer à 100 

for j= L:N1(1)
    
Y(j-L+1)=Q(j);

for k=1:na
Phi((j-L+1),k) =[-Q(j-k)];
end    
    

for k=1:nb
Phi((j-L+1),k+na) =[P(j-k+1-nk)];
end

end


%/////// Conditions intiales  //////////////////////////////
%/// Estimation initiale de la matrice Pm ///////////////////
Pm = inv((Phi')*Phi);
%Pm
%//// Estimation des parametres d'ajustement/////////////////////////
Teta = Pm*(Phi')*Y';
%Teta

%/////////////////////////////////////////////////////////////////////////////
%***************************************************************************** 

%*****************************************************************************
%/////////////////////////////////////////////////////////////////////////////
%/// Algorithme de mise a jour en incorporant la nouvelle information ////////
%/////////////////////////////////////////////////////////////////////////////

function [yh,A]=AlgorMiseAJour(N2,P,Q,na,nb,nk,Pm,Teta,Lamda)

I=eye(na+nb,na+nb);

g=0; %// Cette variable est utilise par la barre de progression

for i= (N2(1)+1):N2(2)   

g=g+1;
%/////// Estimation de Phi //////////////////

for k=1:na        
Ph(k) =[-Q(i-k)]  ;  
end

for k=1:nb
Ph(k+na) =[P(i-k+1-nk)];
end
 
%Ph
%/// Estimation de la matrice des gain K
K=[Pm*(Ph')]/[Lamda+((Ph')')*Pm*(Ph')];
%/// Estimation de la matrice pm
Pm=[I-K*((Ph')')]*(Pm/Lamda);     
%/// Prevision
yh(i)=((Ph')')*Teta;  %// Prevision a priori
%/// Estimation des parametres d'ajustement (Teta) du modele ARX
Teta= Teta + K*[Q(i)-((Ph')')*Teta];
%yh(i)=((Ph')')*Teta;      %// Prevision a posteriori
waitbar(g/((N2(2)-N2(1))))

%// Estimation du modele discret en terme de fonction de transfert H(z)
A(:,i)=Teta;  %// Enregistrement des parametres Teta ds la matrice A
   
end

%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////////
function Dessiner(yh,Q,N1,na,nb,A);


%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees et des previsions //////////////
%//////////////////////////////////////////////////////////////////////////
figure,subplot(2,2,1), plot((N1(1)+1: N1(2)),Q(N1(1)+1: N1(2)));
hold on, plot((N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)),'r');
title('Comparaison entre le débit mesuré et le débit prédit');
legend('Débit mesuré','Débit prédit ARX (MMCRP)') ;
xlabel('Temps en pas d1 heure');
ylabel('Débit (mm)');, grid on

%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des erreurs de previsions //////////////
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))'] ;
subplot(2,2,3), plot((N1(1)+1: N1(2)),100*e/max(Q(N1(1)+1: N1(2))))
ylabel('Erreur de prédiction en %'), grid on
xlabel('Temps en pas d1 jour')
legend('Erreur de prédiction en %') ;

%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees en fonction des previsions /////
%//////////////////////////////////////////////////////////////////////////
subplot(2,2,2),plot(Q(N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)),'o');
title('Relation entre le débit mesuré et le débit prédit');
legend('nuage de points') ;
xlabel('Débit mesuré [mm]');
ylabel('Débit prédit [mm]');

%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des parametres adaptatifs du modele  //////////////
%//////////////////////////////////////////////////////////////////////////
for i= 1:na+nb  %// Balayage de toute les lignes de la matrices A
    m=A(i,:);   %// Chargement du contenu des ligne
     hold on,subplot(2,2,4),plot((N1(1)+1: N1(2)),m((N1(1)+1):N1(2)));
    title('Paramètres identifiés du modèle ARX par MMCRP');
xlabel('Temps: Pas de Temps Journalier');
ylabel('Paramètres du modèle ARX');
end


%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
function Critere_Validation(Q,yh,N1,na,nb)
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))]' ;  %// Calcul de l'erreur de prevision
%/////////////////////////////////////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%////  Autocorrélation de l'erreur de prédiction (MMCRP)  /////////////////
N=length(e);
[R,s]=xcorr(e,'biased');
mx=max(R);
figure,plot(s(N-2:end),R(N-2:end)/mx), grid on
title('Autocorrélation de Erreur de prédiction');
xlabel('Temps: Année 1994 à 1996 avec un pas de temps journalier');
ylabel('Autocorrélation de Erreur de prévision');
sprintf('FONCTION D AUTOCORRELATION: R(0), RN(0)...RN(10)') 
sprintf(' R(0)') 
R(N)
sprintf('RN(0)...RN(10)') 
(R(N:N+10)/mx )' %//
%/////////////  Test de blancheur de l'erreur de prédiction  /////////////
a=1.96/sqrt(N);
f=ones(N);
hold on, plot(a*f(1,:),'r');
hold on, plot(-a*f(1,:),'r');
legend('Autocorrélation','Intervalle de confiance 95%: Borne Sup','Intervalle de confiance 95%: Borne Inf' ) ;
%//////////////////////////////////////////////////////////////////
%/// Fonction d'erreur ///////////////////////////////////////////
%//////////////////////////////////////////////////////////////////
%/// Le biais  //////////////////////////////////////////////////
m1=mean(Q(N1(1)+1: N1(2)));
m2=mean(yh(N1(1)+1:N1(2)));
Biais=100*[(m2-m1)/m1]
%///////////////////////////////////////////////////////////
%/// Coefficient de Nash-Suttcliffe  ///////////////////////
E1=(e'*e);
E2=(Q(N1(1)+1: N1(2))-m1);
E2=E2'*E2;
Nash= 1-(E1/E2) 
%//////////////////////////////////////////////////////////////////////////
%// Coefficient de correlation 
%/////////////////////////////////////////////////////////////////////////
r=corrcoef(Q(N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)));
r=r(1,2)
%/////////////////////////////////////////////////////////////////////////
%*****************************************************************************
%///////////////////////////////////////////////////////////
%/// AIC: Akaike Information Criterion  ////////////////////
E=sum(e);
E=abs(E);
AIC=log10(E)+2*((na+nb)/N)
%/// FPE: Final Prediction Error  ////////////////////
FPE=E1*((N+na+nb-1)/(N-na-nb-1))
%/// Variance des erreurs  ////////////////////
VarE=E1/N
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************

%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************

%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%// Estimation du modele discret en terme de fonction de transfert Q(z)/P(z)=H(z)
%//////////////////////////////////////////////////////////////////////////

function Modele(N1,A,na,nb,Q,P,nk)


for i= 1:na+nb  %// Balayage de toute les lignes de la matrices A
    m=A(i,:);   %// Chargement du contenu des ligne
    m=A(i,(N1(1)+1):N1(2));  
    g(i)=mean(m);  %// Valeur moyenne des parametres
    v(i)=std(m,1); %// Ecart-type (Voir help std)
    
end
g,v

den=[1,-g(1:na)];
num1=[g(na+1:na+nb)];

if  (na>0)&(nb>0)
%//////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////
  sys1 = tf(num1,den,1)
Z1=zpk(sys1)
figure,pzmap(sys1)
pole(sys1)
zero(sys1)
%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
sysc = d2c(sys1,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys1)
figure,step(H)
figure, impulse(H)
figure,pzmap(H)
  
elseif (na==0)&(nb>0)
%//////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////    
sys = tf(num1,1,1)
Z=zpk(sys)
figure,pzmap(sys)
pole(sys)
zero(sys)
%//////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////
sysc = d2c(sys,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys)
end
