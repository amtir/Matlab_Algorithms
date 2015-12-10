

%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%///// Methode des moindres carres  ///////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%// Considerons le modele ARX avec une variable explicative U=P(pr�cipitation)
%// la sortie est le d�bit Q, d'ordre ARX(na,nb) //////////////////////////
%// A(z)Y(k)=B(Z)U(k-nk)+e(k)//////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%// On doit estimer les coefficients ai i=1,2,...,na et bj, j=1,2,...,nb//
%// Donc en tout, na+nb parametres a estimer(inconues) --> na+nb equations minimum //
%//////////////////////////////////////////////////////////////////////////
%// Periode de calibration: 
%// N1=[2300,2500];              
%// Periode de prevision: 
%// N2=[3000,3500]
%// Ordre du systeme: na=3, nb=3, nk=0;
%//////////////////////////////////////////////////////////////////////////
% REFERENCE: R Longchamp, Commande Num�rique de Syst�mes Dynamiques
%                         Presse Polytechniques et universitaires Romandes
%**************************************************************************


function   MMC


%**************************************************************************
%*****************************************************************************
%///// Programme principale //////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%////////  1 Sous fonction : Charger  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[P,Q]=Charger;

%//////////////////////////////////////////////////////////////////////////
%////////  2 Sous fonction : Partie interactive  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
sprintf('Ce programme permet d estimer le vecteur param�tre (MMC) et de simuler le d�bit') 
reply = input('Voulez-Vous continuer? oui=1, non=0 :');
if isempty(reply)
    quit;
elseif reply
N1 = input('Choisir la p�riode de crue pour le calibrage des param�tres du mod�le ARX? : N1=[N1(1),N1(2)]=');
N2 = input('Choisir la p�riode de crue pour la simulation? : N2=[N2(1),N2(2)]=');
na = input('Choisir l ordre pour AR, na=');
nb = input('Choisir l ordre pour la partie explicative, nb=');
nk = input('Choisir le retard, si vous ne savez pas, mettre nk=0, nk=');
H = waitbar(0,'SVP, Attendre la fin des calculs...');

%//////////////////////////////////////////////////////////////////////////
%////////  3 Sous fonction : EstParm  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[Teta]=EstParm(N1,na,nb,Q,P,nk);
abs(Teta(2:end)*1000)
%//////////////////////////////////////////////////////////////////////////
%////////  4 Sous fonction : Simulation  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[yh]=Simulation(N2,P,Q,na,nb,nk,Teta);
            
            
%//////////////////////////////////////////////////////////////////////////
%////////  5 Sous fonction : Dessiner  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
Dessiner(yh,Q,P,N2);
%//////////////////////////////////////////////////////////////////////////
%////////  6 Sous fonction : Critere_Validation  //////////////////////////
%//////////////////////////////////////////////////////////////////////////
Critere_Validation(Q,yh,N2,na,nb);
%//////////////////////////////////////////////////////////////////////////

Modele(Teta,na,nb)

close(H)
end

%*****************************************************************************
%/////////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////////

function [Teta]=EstParm(N1,na,nb,Q,P,nk)


for j=N1(1):N1(2)
    
    Y(j-N1(1)+1)=Q(j);

for k=1:na
Phi((j-N1(1)+1),k) =[-Q(j-k)];
end    
    

for k=1:nb
Phi((j-N1(1)+1),k+na) =[P(j-k+1-nk)];
end



end

%/////////////////////////////////////////////////////////////////////////////
%/// Estimation des parametres ai et bi 
Pm = inv((Phi')*Phi);
%//// Estimation des parametres d'ajustement/////////////////////////
Teta = Pm*(Phi')*(Y')
%/////////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////

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
xlabel('Temps en pas d1 jour: Janvier 1992 � D�cembre 1998')
%/////////////////////////////////////////////////////////////////////////////
%**************************************************************************
%**************************************************************************



%*****************************************************************************
%/////////////////////////////////////////////////////////////////////////////
%/////////////////  Plot des debits mesurees et des previsions //////////////
%/////////////////////////////////////////////////////////////////////////////
function Dessiner(yh,Q,P,N2);

%/////// Plot des Pr�cipitations  //////////////
figure,subplot(4,1,1),plot((N2(1)+1:N2(2)),P(N2(1)+1: N2(2))); 
%title('Precipitation Equivalente ');
ylabel('Precipitation Equivalente'), grid on
xlabel('Temps en pas d1 jour')
legend('Precipitation Equivalente') ;

%/////// Plot des Pr�cipitations /////////////
subplot(4,1,2),plot((N2(1)+1:N2(2)),Q(N2(1)+1: N2(2)),'k');
hold on, plot((N2(1)+1:N2(2)),yh(N2(1)+1:N2(2))','r');
legend('D�bit mesur�','D�bit pr�dit ARX (MMC)') ;
xlabel('Temps: Pas de temps journalier');
ylabel('D�bit');

%/////// Calcul et plot des erreurs  //////////////
e(N2(1)+1: N2(2))=Q(N2(1)+1: N2(2))-(yh(N2(1)+1: N2(2)))';
subplot(4,1,3),plot((N2(1)+1: N2(2)),100*(e(N2(1)+1: N2(2)))/max(Q(N2(1)+1: N2(2)))); 
%title('Erreur de pr�diction en %');
ylabel('Erreur de pr�diction en %'), grid on
xlabel('Temps en pas d1 jour')
legend('Erreur de pr�diction en %') ;


%//////////////////////////////////////////////////////////////////////////
%// Plot du debit mesure en fonction du debit simule
%//////////////////////////////////////////////////////////////////////////
subplot(4,1,4), plot(yh(N2(1)+1:N2(2)),Q(N2(1)+1: N2(2)),'o');
legend('nuage de points') ;
xlabel('D�bit mesur� ');
ylabel('D�bit pr�dit ');


%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************

function [yh]=Simulation(N2,P,Q,na,nb,nk,Teta)
      
%///////////////////////////////////////////////  
for i=N2(1)+1:N2(2)  
    
%/////// Estimation de Phi //////////////////
%/// Vecteur de mesure ou d'observation /////
for k=1:na        
Ph(k) =[-Q(i-k)]  ;  
end

for k=1:nb
Ph(k+na) =[P(i-k+1-nk)];
end
 
%/// Estimation des erreurs de pr�visions du modele ARX
yh(i)=Ph*Teta;   %// Pr�vision

end
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************



function Critere_Validation(Q,yh,N1,na,nb)

%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))]' ;  %// Calcul de l'erreur de prevision
%/////////////////////////////////////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%////  Autocorr�lation de l'erreur de pr�diction (MMCRP)  /////////////////
N=length(e);
[R,s]=xcorr(e,'biased');
mx=max(R);
figure,plot(s(N-2:end),R(N-2:end)/mx), grid on
title('Autocorr�lation de Erreur de pr�diction');
xlabel('Temps: Ann�e 1994 � 1996 avec un pas de temps journalier');
ylabel('Autocorr�lation de Erreur de pr�vision');
sprintf('FONCTION D AUTOCORRELATION: R(0), RN(0)...RN(10)') 
sprintf(' R(0)') 
R(N)
sprintf('RN(0)...RN(10)') 
(R(N:N+10)/mx )' %//
%/////////////  Test de blancheur de l'erreur de pr�diction  /////////////
a=1.96/sqrt(N);
f=ones(N);
hold on, plot(a*f(1,:),'r');
hold on, plot(-a*f(1,:),'r');
legend('Autocorr�lation','Intervalle de confiance 95%: Borne Sup','Intervalle de confiance 95%: Borne Inf' ) ;
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
%//////////////////////////////////////////////////////////////////////////
%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////

function Modele(Teta,na,nb)


den=[1,-Teta(1:na)'];
num1=[Teta(na+1:na+nb)'];

if  (na>0)&(nb>0)
%//////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////
  sys1 = tf(num1,den,1)
  pole(sys1)
Z1=zpk(sys1)
%figure,pzmap(sys1)
%pole(sys1)
%zero(sys1)
%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
sysc = d2c(sys1,'tustin')
H = zpk(sysc)
pole(sysc)
[z,p,k] = zpkdata(sys1)
%figure,step(H)
%figure, impulse(H)
%figure,pzmap(H)
  
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
sysc = d2c(sys) %,'tustin')
H = zpk(sysc)
[z,p,k] = zpkdata(sys)
end


%//////////////////////////////////////////////////////////////////////////