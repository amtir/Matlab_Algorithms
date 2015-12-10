



%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%///// FILTRE NUMERIQUE NON-STATIONAIRE  //////////////////////////////////
%///////////////// FILTRE A MEMOIRE LIMITEE ////////////////////////////////
%//////Methode des moindres carres ponderes recurrentts////////////////////
%// Considerons le modele ARMAX avec une variable explicative U=P(précipitation)
%// la sortie est le débit Q, d'ordre ARMAX(na,nb,nc,nk), par exemple na=10 
%//  nb=5, nc=5, nk=0//////////////////////////////////////////////////////
%// A(z)Y(k)=B(Z)U(k)+e(k) /////////////////////////////////////////////////
%///// Periode de prevision N1=[N1(1), N1(2)] // N1(1)> 350 ////////////////
%//////////// nk est le retard de la variable explictive////////////////////
%/////////// 0.95<Lamda<0.99 est le facteur d'oubli ////////////////////////
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%// Modele ARMAX, ARX, AR, ARMA, MA à paramètres Adaptatifs ///////////////
%//////////////////////////////////////////////////////////////////////////

function MMCRPE_Matlab

%N1=[800,1500]; na=5, nb=5,nk=0, Lamda=0.999;
%**************************************************************************
%////////// Programme principale ///////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%////////  1 Sous fonction : Charger  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[P,Q]=Charger;

%//////////////////////////////////////////////////////////////////////////
%////////  2 Sous fonction : Partie interactive  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
sprintf('Ce programme permet d estimer le vecteur paramètre et de simuler le débit') 
reply = input('Voulez-Vous continuer? oui=1, non=0 :');
if isempty(reply)
    quit;
elseif reply
N1 = input('Choisir la période de crue pour la simulation? : N1=[N1(1),N1(2)]=');
na = input('Choisir l ordre pour AR, na=');
nb = input('Choisir l ordre pour la partie explicative, nb=');
nc = input('Choisir l ordre pour la partie MA, nc=');
nk = input('Choisir le retard, si vous ne savez pas, mettre nk=0, nk=');
Lamda = input('Choisir le coefficient de pondération entre [0.98,0.99], Lamda=');
H = waitbar(0,'SVP, Attendre la fin des calculs...');
end

%[Q]=Methode_Lissage(Q,20,1,0.9);   %// Filtrage du debit 


if (na>0&nb==0&nc==0)  %// Modele AR: Ay(t)=e(t)
     sprintf('Moddele AR') 
elseif (na==0 & nb>0&nc==0)  %// y(t)=Bu(t)+e(t) 
     sprintf('Moddele X')    
elseif (na>0&nb>0&nc==0)  %// Modele ARMAX,ARX, AR, MAX, X
     sprintf('Moddele ARX') 
elseif (na>0&nb==0&nc>0)  %// Modele ARMAX,ARX, AR, MAX, X
     sprintf('Moddele ARMA') 
elseif (na==0&nb==0&nc==0)
    sprintf('Erreur!!! Choisir l Ordre: na, nb, nc >0.')   
elseif (na>0&nb>0&nc>0)
     sprintf('Moddele ARMAX')  
end

%if (na>0&nb>0&nc>0)  %// Modele ARMAX,ARX, AR, MAX, X

%//////////////////////////////////////////////////////////////////////////
%////////  3 Sous fonction : ConditionInitiale0  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[Pm0,Teta0]=ConditionInitiale0(na,nb,nk,Q,P,N1);
%//////////////////////////////////////////////////////////////////////////
%////////  4 Sous fonction : Erreur  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[e]=Erreur(N1,P,Q,na,nb,nk,Pm0,Teta0,Lamda);  %/// Utilisation de ces erreurs pour amorcer l'algorithme 
%//////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%////////  5 Sous fonction : ConditionInitiale  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[Pm,Teta]=ConditionInitiale(na,nb,nc,nk,N1,Q,P,e);
%//////////////////////////////////////////////////////////////////////////
%////////  6 Sous fonction : AlgorithmeMiseAjour  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%[yh,A]=AlgorMiseAJour(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda);
[yh,A]=AlgorMiseAJour_Matlab(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda);
%//////////////////////////////////////////////////////////////////////////
%////////  7 Sous fonction : Dessiner  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
Dessiner(yh,Q,N1,na,nb,nc,A);
%//////////////////////////////////////////////////////////////////////////
%////////  8 Sous fonction : Critere_Validation  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
Critere_Validation(Q,yh,N1,na,nb,nc);
%//////////////////////////////////////////////////////////////////////////
%////////  9 Sous fonction : Modele  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%Modele(N1,A,na,nb,nc,Q,P,nk);
close(H)

%elseif (na>0|nb>0&nc==0)
 %   e=0;
%[Pm,Teta]=ConditionInitiale(na,nb,nc,nk,N1,Q,P,e);
%[yh,A]=AlgorMiseAJour(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda);
%Dessiner(yh,Q,N1,na,nb,nc,A);
%Critere_Validation(Q,yh,N1,na,nb,nc);
%Modele(N1,A,na,nb,nc,Q,P,nk);
%close(H)

%elseif (na==0&nb==0&nc>0)  %// Modele MA  
%sprintf('Moddele MA') 
%e=Q/20;   %// Choix arbitraire des erreurs
%[Pm,Teta]=ConditionInitiale(na,nb,nc,nk,N1,Q,P,e);
%[yh,A]=AlgorMiseAJour(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda);
%N1=[N1(1)+20,N1(2)];
%Dessiner(yh,Q,N1,na,nb,nc,A);
%Critere_Validation(Q,yh,N1,na,nb,nc);
%Modele(N1,A,na,nb,nc,Q,P,nk);
%close(H)
%end
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*************************************************************************
%//////////////////////////////////////////////////////////////////////////
%Commencer par charger les données: moyenne journalière du bassin /////////
function [P,Q]=Charger


%// Charger les Pluie equivalente de la partie non glaciere (fonte + pluie)
load('Peq_glacier_gletsch_journalier_92_2001.txt');  
Peq_ng=Peq_glacier_gletsch_journalier_92_2001(2:end);

%// Charger les Pluie equivalente de la partie glaciere (pluie +fonte neige)
load('Peq_non_glacier_gletsch_journalier_92_2001.txt ');  
Peq_g=Peq_non_glacier_gletsch_journalier_92_2001(2:end);

%// Charger les fonte de du glacier lorsqu'il n'y a pas de neige
load('Mgl_glacier_gletsch_journalier_92_2001.txt');  
Pmgl=Mgl_glacier_gletsch_journalier_92_2001(2:end);

P=(Peq_ng+Peq_g+Pmgl);   %// la pluie éauivalente totale est la somme des 3 contributions

load('debits_gletsh_journalier_92_2001.txt'); %Débits à l'exutoire des deux bassins (Gletsh/Alptal)
Q=debits_gletsh_journalier_92_2001(:,1);  %Débits du Bassin (1)Alptal, (2)Gletsch

figure,subplot(2,1,1), plot(Q), 
ylabel('Debits à l exutoire du Bassin Gletsch (mm)'), grid on
subplot(2,1,2), plot(P)
ylabel('Precipitation Effective (mm)'), grid on
xlabel('Temps en pas d1 jour: De Janvier 1992 au 31 décembre 2001')
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


function [Pm0,Teta0]=ConditionInitiale0(na,nb,nk,Q,P,N1)


L=100;

for j=L:N1(1)
    Y(j-L+1)=Q(j);
    
for i=1:na
    Phi(j-L+1,i)=[-Q(j-i)];
end 

for i=1:nb
    Phi(j-L+1,i+na)=[P(j-i+1-nk)];
end

end

%/////// Conditions intiales  //////////////////////////////
%/// Estimation initiale de la matrice Pm ///////////////////
Pm0 = inv((Phi')*Phi);
Teta0 = Pm0*(Phi')*Y';
%/////////////////////////////////////////////////////////////////////////////
%***************************************************************************** 
%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
%//// Conditions initiales: !!! Informations extraite de l'historique//////
%//////////////////////////////////////////////////////////////////////////

function [Pm,Teta]=ConditionInitiale(na,nb,nc,nk,N1,Q,P,e);


L=100;

for j=L:N1(1)
    Y(j-L+1)=Q(j);
    
for i=1:na
    Phi(j-L+1,i)=[-Q(j-i)];
end 

for i=1:nb
    Phi(j-L+1,i+na)=[P(j-i-nk+1)];
end

for i=1:nc
    Phi(j-L+1,i+na+nb)=e(j-i);
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

function [e]=Erreur(N1,P,Q,na,nb,nk,Pm,Teta,Lamda)

I=eye(na+nb,na+nb);

for i= 80:N1(1)

%/////// Estimation de Phi //////////////////

for k=1:na        
Ph(k) =[-Q(i-k)]  ;  
end

for k=1:nb
Ph(k+na) =[P(i-k-nk+1)];
end


%/// Estimation de la matrice des gain K
%K=[Pm*(Ph')]/[Lamda+((Ph')')*Pm*(Ph')];
    %K=[Pm*(Ph')]/[(1/Lamda)+((Ph')')*Pm*(Ph')];
%/// Estimation de la matrice pm
    %Pm=[I-K*((Ph')')]*(Pm/Lamda);    
%Pm=[I-K*((Ph')')]*Pm;  
%/// Estimation des parametres d'ajustement (Teta) du modele ARMAX, ARX, AR
e(i)=[Q(i)-((Ph')')*Teta]; %/[Lamda+((Ph')')*Pm*(Ph')];  %// Calcul de l'erreur
%Teta= Teta + K*[Q(i)-((Ph')')*Teta];

end
%plot(e)
%*****************************************************************************
%/////////////////////////////////////////////////////////////////////////////
%/// Algorithme de mise a jour en incorporant la nouvelle information ////////
%/////////////////////////////////////////////////////////////////////////////


function [yh,A]=AlgorMiseAJour_Matlab(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda)

nn = [na nb nc nk];
    %nf=nb; nn=[na,nf,nk];
z=[Q(N1(1):N1(2)),P(N1(1):N1(2))];
[thm,yh,P] = rarmax(z,nn,'ff',Lamda,Teta,Pm);
    %[thm,yh,P] = roe(z,nn,'ff',Lamda,Teta,Pm);
figure,plot(thm);
figure,plot((N1(1):N1(2)),yh,'r');
hold on, plot((N1(1):N1(2)),Q(N1(1):N1(2)));
A=thm';
yh(N1(1):N1(2))=yh;
%yh=yh';

function [yh,A]=AlgorMiseAJour(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda)

I=eye(na+nb+nc,na+nb+nc);
g=0; %// Cette variable est utilise pour montrer la 


for i= (N1(1)+1):N1(2)   
g=g+1;
%/////// Estimation de Phi //////////////////

for k=1:na        
Ph(k) =[-Q(i-k)]  ;  
end

for k=1:nb
Ph(k+na) =[P(i-k-nk+1)];
end
 
for k=1:nc
Ph(k+na+nb) =[e(i-k)];
end
%Ph
%/// Estimation de la matrice des gain K
K=[Pm*(Ph')]/[Lamda+((Ph')')*Pm*(Ph')];
%/// Estimation de la matrice pm
Pm=[I-K*((Ph')')]*(Pm/Lamda);     
%/// Estimation des parametres d'ajustement (Teta) du modele ARX
e(i)=[Q(i)-((Ph')')*Teta]/[Lamda+((Ph')')*Pm*(Ph')]; %//  Prendre l'erreur à postériori --> meilleur convergence!!!
%// Ioan Doré Landau, identification et commande, Hermes
yh(i)=((Ph')')*Teta;
Teta= Teta + K*[Q(i)-((Ph')')*Teta];
    %yh(i)=((Ph')')*Teta;
%e
waitbar(g/((N1(2)-N1(1))))
%// Estimation du modele discret en terme de fonction de transfert H(z)
A(:,i)=Teta;  %// Enregistrement des parametres Teta ds la matrice A

end



%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*****************************************************************************
%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees et des previsions //////////////
%/////////////////////////////////////////////////////////////////////////////
function Dessiner(yh,Q,N1,na,nb,nc,A)


%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees et des previsions //////////////
%//////////////////////////////////////////////////////////////////////////
figure,subplot(2,2,1), plot((N1(1)+1: N1(2)),Q(N1(1)+1: N1(2)));
hold on, plot((N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)),'r');
title('Comparaison entre le débit mesuré et le débit prédit');
legend('Débit mesuré','Débit prédit') ;
xlabel('Temps en pas d1 jour');
ylabel('Débit m3/s');, grid on

%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des erreurs de previsions //////////////
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))] ;
subplot(2,2,3), plot((N1(1)+1: N1(2)),100*e/max(Q(N1(1)+1: N1(2))));
ylabel('Erreur de prédiction en %'), grid on;
xlabel('Temps en pas d1 jour');
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
%for i= 1:na+nb+nc %// Balayage de toute les lignes de la matrices A
 %   m=A(i,:);   %// Chargement du contenu des ligne
  %   hold on,subplot(2,2,4),plot((N1(1)+1: N1(2)),m((N1(1)+1):N1(2)));
   % title('Paramètres identifiés du modèle ARX par MMCRP');
%xlabel('Temps: Pas de Temps Journalier');
%ylabel('Paramètres du modèle ARX');
%end

function Critere_Validation(Q,yh,N1,na,nb,nc)
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))] ;  %// Calcul de l'erreur de prevision
%/////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%////  Autocorrélation de l'erreur de prédiction (MMCRP)  /////////////////
N=length(e)
[R,s]=xcorr(e,'biased');
mx=max(R);
figure,plot(s(N-2:end),R(N-2:end)/mx), grid on
title('Autocorrélation de Erreur de prédiction');
xlabel('Temps: Année 1994 à 1996 avec un pas de temps journalier');
ylabel('Autocorrélation de Erreur de prévision');
%/////////////  Test de blancheur de l'erreur de prédiction  /////////////
a=1.96/sqrt(N);
f=ones(N);
hold on, plot(a*f(1,:),'r');
hold on, plot(-a*f(1,:),'r');
legend('Autocorrélation','Intervalle de confiance 95%: Borne Sup','Intervalle de confiance 95%: Borne Inf' ) ;
%//////////////////////////////////////////////////////////////////
%/// Fonctio0n d'erreur ///////////////////////////////////////////
%//////////////////////////////////////////////////////////////////
%/// Le biais  //////////////////////////////////////////////////
m1=mean(Q(N1(1)+1: N1(2)));
m2=mean(yh(N1(1)+1:N1(2)));
Biais=100*[(m2-m1)/m1]
%///////////////////////////////////////////////////////////
%/// Coefficient de Nash-Suttcliffe  ///////////////////////
 E1=(e'*e);
    E2=(Q(N1(1)+2: N1(2))-m1);
E2=E2'*E2;
Nash= 1-(E1/E2) 
%/////////////////////////////////////////////////////////////////////////
%// Coefficient de correlation 
%/////////////////////////////////////////////////////////////////////////
r=corrcoef(Q(N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)));
r=r(1,2)
%//////////////////////////////////////////////////////////////////////////
%*****************************************************************************
%///////////////////////////////////////////////////////////
%/// AIC: Akaike Information Criterion  ////////////////////
E=sum(e);
E=abs(E);
AIC=log10(E)+2*((na+nb+nc)/N)
%/// FPE: Final Prediction Error  ////////////////////
FPE=E1*((N+na+nb+nc-1)/(N-na-nb-nc-1))
%/// Variance des erreurs  ////////////////////
VarE=E1/N
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************
%**************************************************************************
%// Estimation du modele discret en terme de fonction de transfert Q(z)/P(z)=H(z)

function Modele(N1,A,na,nb,nc,Q,P,nk)

%figure, grid on;
%title('Paramètres identifiés du modèle ARX par MMCRP');
%xlabel('Temps: Année 1994 à 1996 avec un pas de temps journalier');
%ylabel('Paramètres du modèle ARX');

for i= 1:na+nb+nc  %// Balayage de toute les lignes de la matrices A
    m=A(i,:);   %// Chargement du contenu des ligne
    m=A(i,(N1(1)+1):N1(2));  
    g(i)=mean(m);  %// Valeur moyenne des parametres
    v(i)=std(m,1); %// Ecart-type (Voir help std)
    
end
g,v

den=[1,-g(1:na)];
num1=[g(na+1:na+nb)];
num2=[g(nb+1:na+nb+nc)];


if (na>0)&(nb>0)&(nc>0)
sys1 = tf(num1,den,1)
Z1=zpk(sys1)
figure,pzmap(sys1)
pole(sys1)
zero(sys1)
%step(Z)
%figure, impulse(Z)
sysc = d2c(sys1,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys1)


sys2 = tf(num2,den,1)
Z2=zpk(sys2)
figure,pzmap(sys2)
pole(sys2)
zero(sys2)
%step(Z)
%figure, impulse(Z)
sysc = d2c(sys2,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys2)

elseif (na>0)&(nb>0)&(nc==0)
  sys1 = tf(num1,den,1)
Z1=zpk(sys1)
figure,pzmap(sys1)
pole(sys1)
zero(sys1)
%step(Z)
%figure, impulse(Z)
sysc = d2c(sys1,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys1)
figure,step(H)
figure, impulse(H)
figure,pzmap(H)
%t=[1:20];
%u=[1:20];
%u(1:5)=0;
%u(6:14)=1;
%u(15:20)=0
%figure,lsim(H,u,t)
  elseif (na>0)&(nb==0)&(nc>0)
sys2 = tf(num2,den,1)
Z2=zpk(sys2)
figure,pzmap(sys2)
pole(sys2)
zero(sys2)
%step(Z)
%figure, impulse(Z)
sysc = d2c(sys2,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys2)

elseif (na==0)&(nb==0)&(nc>0)
sys = tf(num2,1,1)
Z=zpk(sys)
figure,pzmap(sys)
pole(sys)
zero(sys)
%step(Z)
%figure, impulse(Z)
sysc = d2c(sys,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys)

elseif (na==0)&(nb>0)&(nc==0)
sys = tf(num1,1,1)
Z=zpk(sys)
figure,pzmap(sys)
pole(sys)
zero(sys)
%step(Z)
%figure, impulse(Z)
sysc = d2c(sys,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys)
end
%m=size(p)
%a=[1;1]
%a=p{1,1}
%A=k/(a(1)-a(2))
%B=k/(a(2)-a(1))
%a1=[a(1)]
%b1=[a(2)]
%Tet1=[b1,B]
%Tet2=[a1,A]
%for i= (N1(1)+1):N1(2)


%for j=1:na-1
 %   Z(j)=[-Q(i-j)]
 %  end

    %for k=1:nb1
    %pr(k)=[P(i-nk-k)]
    %end  
    
    %N=[Z,pr];
%yh2(i+1) = N*Tet1' ;
%yh3(i+1) = N*Tet2' ;
%end

%figure(1),hold on, plot(yh2(N1(1)+2:N1(2)),'g');
%hold on, plot(yh3(N1(1)+2:N1(2)),'r');