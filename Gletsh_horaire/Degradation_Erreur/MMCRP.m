



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


function [e,yh,Q,u,U,St]=MMCRP(N1,na,nb,nk,Lamda)


%**************************************************************************
%////////// Programme principale ///////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%////////  1 Sous fonction : Charger  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[P,Q]=Charger;


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
[e,u,U,St]=Critere_Validation(Q,yh,N1,na,nb);
%//////////////////////////////////////////////////////////////////////////
%////////  7 Sous fonction : Modele  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%Modele(N1,A,na,nb,Q,P,nk);


%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*************************************************************************
%//////////////////////////////////////////////////////////////////////////
%Commencer par charger les données: moyenne journalière du bassin /////////
function [P,Q]=Charger

load('Gletch.txt');  %Charger les inputs

Peq_non_glacier_gletsch=Gletch(1:5136)  ;
Peq_glacier_gletsch=Gletch(5137:10272) ;
Mgl_glacier_gletsch=Gletch(10273:15408)  ;

P=Peq_non_glacier_gletsch+Peq_glacier_gletsch+Mgl_glacier_gletsch ;
load('aaDEBIT_horaire2.txt');
%Charger les donnees horaires du bassin Gletsch
D=aaDEBIT_horaire2(1:5136);
%charger les donnees horaires de bassin Aptal
D2=aaDEBIT_horaire2(10273:19032);
%Charger les inputs horaires du bassin Gletsch (pluies equivalentes)
%load('Gletch.txt');
%P=Gletch(10273:15408);
Q=D;

%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************

%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
%//// Amorçage de l'algorithme à partir de N2  ////////////////////////////
%//////////////////////////////////////////////////////////////////////////

function [Pm,Teta]=ConditionInitiale(na,nb,nk,N1,Q,P);


%L=max(na,(nb+nk))+1;  %// Utiliser le maximum d'historique
%L=N1(1)-400;          %// Prendre une année (400jour) d'historique
L=100;                 %// Commencer à 100 

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

yh(N2(1))=Q(N2(1));

I=eye(na+nb,na+nb);
yh(N2(1))=Q(N2(1));

for i= (N2(1)+1):N2(2)   


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
Q(i)=((Ph')')*Teta;   %// Estimation du débit mesuré par la valeur prédite 
%/// Estimation des parametres d'ajustement (Teta) du modele ARX
Teta= Teta + K*[Q(i)-((Ph')')*Teta];
%yh(i)=((Ph')')*Teta;      %// Prevision a posteriori



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
legend('Débit mesuré','Débit prédit ARX (MMCRPE)') ;
xlabel('Temps en pas d1 jour');
ylabel('Débit m3/s');, grid on

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
xlabel('Débit mesuré m3/s');
ylabel('Débit prédit m3/s');

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
function [epct,u,U,St]=Critere_Validation(Q,yh,N1,na,nb)
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1): N1(2))]- [yh(N1(1):N1(2))]' ;  %// Calcul de l'erreur de prevision
%/////////////////////////////////////////////////////////////////////////
%epct=100*(e./mean(Q(N1(1)+1: N1(2))));  %Erreur de prédiction en %
epct=e./10;
fin=length(yh);
St(N1(1))=0;
for i=N1(1):fin
    Y=yh(N1(1):i);
St(i)=std((Y./10),1);
end

u=100*(e./Q(N1(1): N1(2)));  %// Erreur de prédiction en %
fin=length(u);
U(1)=u(1);
for i=2:fin
    U(i)=U(i-1)+u(i);
end
U=U';

%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%////  Autocorrélation de l'erreur de prédiction (MMCRP)  /////////////////
N=length(e);

%//////////////////////////////////////////////////////////////////
%/// Fonction d'erreur ///////////////////////////////////////////
%//////////////////////////////////////////////////////////////////
%/// Le biais  //////////////////////////////////////////////////
m1=mean(Q(N1(1)+1: N1(2)));
m2=mean(yh(N1(1)+1:N1(2)));
Biais=100*[(m2-m1)/m1];
%///////////////////////////////////////////////////////////
%/// Coefficient de Nash-Suttcliffe  ///////////////////////
E1=(e'*e);
E2=(Q(N1(1)+1: N1(2))-m1);
E2=E2'*E2;
Nash= 1-(E1/E2) ;
%//////////////////////////////////////////////////////////////////////////
%// Coefficient de correlation 
%/////////////////////////////////////////////////////////////////////////
r=corrcoef(Q(N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)));
r=r(1,2);
%/////////////////////////////////////////////////////////////////////////
%*****************************************************************************
%///////////////////////////////////////////////////////////
%/// AIC: Akaike Information Criterion  ////////////////////
E=sum(e);
E=abs(E);
AIC=log10(E)+2*((na+nb)/N);
%/// FPE: Final Prediction Error  ////////////////////
FPE=E1*((N+na+nb-1)/(N-na-nb-1));
%/// Variance des erreurs  ////////////////////
VarE=E1/N;
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
