


%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%///// A.A.P : Algorithme d'adaptation param�trique: M_V_R  ///////////////
%//////////////////////////////////////////////////////////////////////////
%///// Maximum de vraissemblance r�cursif   ///////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%///// FILTRE NUMERIQUE NON-STATIONAIRE, param�tres adaptatifs ////////////
%///// FILTRE A MEMOIRE LIMITEE ///////////////////////////////////////////
%///// Estimation des param�tres du mod�le ARMAX: A(z)Y(k)=B(Z)U(k)+C(z)e(k)
%///// une variable explicative U=P(pr�cipitation)  ///////////////////////
%///// la sortie est le d�bit Q, d'ordre ARMAX(na,nb,nc,nk)////////////////
%///// Periode de prevision N1=[N1(1),N1(2)]///////////////////////////////
%///// na, nb et nc l'ordre du model //////////////////////////////////////
%///// nk est le retard de la variable explicative  ///////////////////////
%///// 0.95<Lamda<0.99 est le facteur d'oubli /////////////////////////////
%///// Exemple na=3, nb=3, nc=3, nk=0 et N1=[2000,3500]  //////////////////
%//////////////////////////////////////////////////////////////////////////
% REFERENCE: Ioan D Landau, Identification des syst�mes
%                           Paris, Herm�s Science publication 2001
%**************************************************************************

function M_V_R

%**************************************************************************
%////////// Programme principale ///////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%////////  1 Sous fonction : Charger  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[P,Q]=Charger;

%//////////////////////////////////////////////////////////////////////////
%////////  2 Sous fonction : Partie interactive  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
sprintf('Ce programme permet d estimer le vecteur param�tre et de simuler le d�bit') 
reply = input('Voulez-Vous continuer? oui=1, non=0 :');
if isempty(reply)
    quit;
elseif reply
N1 = input('Choisir la p�riode de crue pour la simulation? : N1=[N1(1),N1(2)]=');
na = input('Choisir l ordre pour AR, na=');
nb = input('Choisir l ordre pour la partie explicative, nb=');
nc = input('Choisir l ordre pour la partie MA, nc=');
nk = input('Choisir le retard, si vous ne savez pas, mettre nk=0, nk=');
Lamda = input('Choisir le coefficient de pond�ration entre [0.98,0.99], Lamda=');
H = waitbar(0,'SVP, Attendre la fin des calculs...');


%[Q]=Methode_Lissage(Q,20,1,0.9);   %// Filtrage du debit 


if (na>0 & nb==0&nc==0)  %// Modele AR: Ay(t)=e(t)
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
%// Estimation des erreurs 
[Pm0,Teta0]=ConditionInitiale0(na,nb,nk,N1,Q,P);
[e]=Erreur(N1,P,Q,na,nb,nk,Pm0,Teta0,Lamda);  %/// Utilisation de ces erreurs pour amorcer l'algorithme 
%//////////////////////////////////////////////
[Pm,Teta]=ConditionInitiale(na,nb,nc,nk,N1,Q,P,e);
[yh,A]=AlgorMiseAJour(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda);
Dessiner(yh,Q,N1,na,nb,nc,A);
Critere_Validation(Q,yh,N1,na,nb,nc)
%Modele(N1,A,na,nb,nc,Q,P,nk);
close(H)
end

%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*************************************************************************
%//////////////////////////////////////////////////////////////////////////
%Commencer par charger les donn�es: moyenne journali�re du bassin /////////
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

P=(Peq_ng+Peq_g+Pmgl);   %// la pluie �auivalente totale est la somme des 3 contributions

load('debits_gletsh_journalier_92_2001.txt'); %D�bits � l'exutoire des deux bassins (Gletsh/Alptal)
Q=debits_gletsh_journalier_92_2001(:,1);  %D�bits du Bassin (1)Alptal, (2)Gletsch



figure,subplot(2,1,1), plot(Q), 
ylabel('Debits � l exutoire du Bassin Gletsch (mm)'), grid on
subplot(2,1,2), plot(P)
ylabel('Precipitation Effective (mm)'), grid on
xlabel('Temps en pas d1 jour: De Janvier 1992 au 31 d�cembre 2001')

%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


function [Pm0,Teta0]=ConditionInitiale0(na,nb,nk,N1,Q,P);

L=100;

for j=L:N1(1)
    Y(j-L+1)=Q(j);
    
for i=1:na
    Phi(j-L+1,i)=[-Q(j-i)];
end 

for i=1:nb
    Phi(j-L+1,i+na)=[P(j-i-nk+1)];
end

end


%/////// Conditions intiales  //////////////////////////////
%/// Estimation initiale de la matrice Pm ///////////////////
Pm0 = inv((Phi')*Phi);
%Pm0
%//// Estimation des parametres d'ajustement/////////////////////////
Teta0 = Pm0*(Phi')*Y';
%Teta0

%/////////////////////////////////////////////////////////////////////////////
%***************************************************************************** 
%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
%//// Conditions initiales: !!! Informations extraite de l'historique//////
%//////////////////////////////////////////////////////////////////////////

function [Pm,Teta]=ConditionInitiale(na,nb,nc,nk,N1,Q,P,e);

i=N1(1);

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

function [yh,A]=AlgorMiseAJour(N1,P,Q,e,na,nb,nc,nk,Pm,Teta,Lamda)

I=eye(na+nb+nc,na+nb+nc);
g=0; %// Cette variable est utilise pour montrer la 

%///////////////////////////////////////
%// Estimation du polynome C
Cf(1:nc+1)=1;
Cf(2:nc+1)=[Teta(na+nb+1:na+nb+nc)];
%// Stabilit� du filtre C 
Cf2=polystab(Cf);
%//////////////////////////////////////
i=N1(1)+1;

for k=1:na     
Q(i-k)=filter(1,Cf2,Q(i-k));  
end

for k=1:nb
P(i-k)=filter(1,Cf2,P(i-k-nk+1));  
end
 
for k=1:nc
e(i-k)=filter(1,Cf2,e(i-k));  
end
%/////////////////////////////////

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


%********************************************************
%// Maximum de vraisemblance r�cursif (M_V_R)
%********************************************************
%// Estimation du polynome C
C(1:nc+1)=1;
C(2:nc+1)=[Teta(na+nb+1:na+nb+nc)];
%// Stabilit� du filtre C 
C2=polystab(C);
%*******************************************************
%//////////////////////////////////////////////////////////////////   
yh(i)=((Ph')')*Teta;  %// Estimation � prioiri
%/// Estimation des parametres d'ajustement (Teta) du modele ARX
e(i)=[Q(i)-((Ph')')*Teta]/[1+((Ph')')*Pm*(Ph')]; %//  Prendre l'erreur � post�riori --> meilleur convergence!!!

Teta=Teta+(Pm*(Ph')*e(i));

%/// Estimation de la matrice pm=F 
Pm=[Pm-((Pm*(Ph')*((Ph')')*Pm)/[Lamda+((Ph')')*Pm*(Ph')])]*(1/Lamda);  
%// Ioan Dor� Landau, identification et commande, Hermes
%///////////////////////////////////////////////////////

A(:,i)=Teta;  %// Enregistrement des parametres Teta ds la matrice A

%// Filtrer le vecteur d'observation avec le filtre 1/C
e(i)=filter(1,C2,e(i));  %// Erreur filtr�e
Q(i)=filter(1,C2,Q(i));  %// Vecteur d'observation filtr�
P(i)=filter(1,C2,P(i));  %// Vecteur d'observation filtr�
waitbar(g/((N1(2)-N1(1))))

end



%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*****************************************************************************
%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees et des previsions //////////////
%/////////////////////////////////////////////////////////////////////////////
function Dessiner(yh,Q,N1,na,nb,nc,A);


%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees et des previsions //////////////
%//////////////////////////////////////////////////////////////////////////
figure,subplot(2,2,1), plot((N1(1)+1: N1(2)),Q(N1(1)+1: N1(2)));
hold on, plot((N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)),'r');
title('Comparaison entre le d�bit mesur� et le d�bit pr�dit');
legend('D�bit mesur�','D�bit pr�dit (MVR)') ;
xlabel('Temps en pas d1 jour');
ylabel('D�bit (mm)');, grid on

%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des erreurs de previsions //////////////
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))'] ;
subplot(2,2,3), plot((N1(1)+1: N1(2)),100*e/max(Q(N1(1)+1: N1(2))))
ylabel('Erreur de pr�diction en %'), grid on
xlabel('Temps en pas d1 jour')
legend('Erreur de pr�diction en %') ;


%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees en fonction des previsions /////
%//////////////////////////////////////////////////////////////////////////
subplot(2,2,2),plot(Q(N1(1)+1: N1(2)),yh(N1(1)+1:N1(2)),'o');
title('Relation entre le d�bit mesur� et le d�bit pr�dit');
legend('nuage de points') ;
xlabel('D�bit mesur� m3/s');
ylabel('D�bit pr�dit m3/s');

%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des parametres adaptatifs du modele  //////////////
%//////////////////////////////////////////////////////////////////////////
for i= 1:na+nb+nc %// Balayage de toute les lignes de la matrices A
    m=A(i,:);   %// Chargement du contenu des ligne
     hold on,subplot(2,2,4),plot((N1(1)+1: N1(2)),m((N1(1)+1):N1(2)));
    title('Param�tres identifi�s du mod�le ARX par MMCRP');
xlabel('Temps: Pas de Temps Journalier');
ylabel('Param�tres du mod�le ARX');
end
%//////////////////////////////////////////////////////////////////////////
function Critere_Validation(Q,yh,N1,na,nb,nc)
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))]' ;  %// Calcul de l'erreur de prevision
%/////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%////  Autocorr�lation de l'erreur de pr�diction (MMCRP)  /////////////////
N=length(e)
[R,s]=xcorr(e,'biased');
mx=max(R);
figure,plot(s(N-2:end),R(N-2:end)/mx), grid on
title('Autocorr�lation de Erreur de pr�diction');
xlabel('Temps: Ann�e 1994 � 1996 avec un pas de temps journalier');
ylabel('Autocorr�lation de Erreur de pr�vision');
%/////////////  Test de blancheur de l'erreur de pr�diction  /////////////
a=1.96/sqrt(N);
f=ones(N);
hold on, plot(a*f(1,:),'r');
hold on, plot(-a*f(1,:),'r');
legend('Autocorr�lation','Intervalle de confiance 95%: Borne Sup','Intervalle de confiance 95%: Borne Inf' ) ;
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
%*****************************************************************************
%**************************************************************************
%// Estimation du modele discret en terme de fonction de transfert Q(z)/P(z)=H(z)

function Modele(N1,A,na,nb,nc,Q,P,nk)

%figure, grid on;
%title('Param�tres identifi�s du mod�le ARX par MMCRP');
%xlabel('Temps: Ann�e 1994 � 1996 avec un pas de temps journalier');
%ylabel('Param�tres du mod�le ARX');

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