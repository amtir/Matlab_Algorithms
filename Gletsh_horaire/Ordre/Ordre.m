
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%///// Programme permettant de determiner l'ordre du modele   /////////////
%///// En analysant les graphes obtenus:                      /////////////
%///// Variance de l'erreur de prediction en fonction de na et nb /////////
%///// na, nb est l'ordre d'un modele ARX /////////////////////////////////
%///// Exemple: N1=[2800,3500]; na=5, nb=5,nk=0      //////////////////////
%**************************************************************************
function Ordre

%**************************************************************************
%////////// Programme principale ///////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////
%////////  1 Sous fonction : Charger  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
[P,Q]=Charger;

%//////////////////////////////////////////////////////////////////////////
%////////  2 Sous fonction : Partie interactive  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
sprintf('Ce programme permet de déterminer l ordre du modèle ARX, na, nb>0') 
reply = input('Voulez-Vous continuer? oui=1, non=0 :');
if isempty(reply)
    quit;
elseif reply
N1 = input('Choisir la période de crue? : N1=[N1(1),N1(2)]=');
na = input('Choisir l ordre max pour AR, na=');
nb = input('Choisir l ordre max pour la partie explicative, nb=');
nk = input('Choisir le retard, si vous ne savez pas, mettre nk=0, nk=');
H = waitbar(0,'SVP, Attendre la fin des calculs...');

%//////////////////////////////////////////////////////////////////////////
%////////  Repetiter des differentes sous fonctions ci-desous  //////////
%/////////        pour different ordre    /////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
for i=1:na   %// different ordre de na
    waitbar(i/na)
    for j=1:nb   %// different ordre de nb
        
        if na==0&nb==0
            sprintf('Erreur!!! Choisir l Ordre: na, nb>0.')   
    
        elseif na>0|nb>0  %// Modele ARX, AR, MAX, X
%//////////////////////////////////////////////////////////////////////////
%////////  3 Sous fonction : ConditionInitiale  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////           
            [Teta]=ConditionInitiale(i,j,nk,N1,Q,P);
%//////////////////////////////////////////////////////////////////////////
%////////  4 Sous fonction : Simulation  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
            [yh]=Simulation(N1,P,Q,i,j,nk,Teta);
%//////////////////////////////////////////////////////////////////////////
%////////  5 Sous fonction : Critere  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
            [VarE,AIC]=Critere(yh,Q,N1,i,j);
        end

Var(i,j)=VarE;
Aic(i,j)=AIC;

    end
    
end
%//////////////////////////////////////////////////////////////////////////
%////////  6 Sous fonction : Dessiner  /////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
Dessiner(Var)
%//////////////////////////////////////////////////////////////////////////
close(H)
end

%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%//////////  Programme des sous fonctions  ////////////////////////////////
%//////////////////////////////////////////////////////////////////////////

function Dessiner(Var)
figure,surf(Var);
legend('Variance de l erreur de prédiction')
xlabel('Ordrede la partie AR: na')
ylabel('Ordre de la partie exogène: nb')
zlabel('Variance de l erreur de prévision')
m1=min(min(Var));   %// Valeur minimale de la variance des erreurs
m2=max(max(Var));   %// Valeur maximale de la variance des erreurs
caxis auto
colorbar('horiz')

%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
figure, subplot(2,1,1), contourf(Var,30)
caxis auto
colorbar('horiz')
xlabel('Ordrede la partie AR: na')
ylabel('Ordre de la partie exogène: nb')
legend('Variance de l erreur de prédiction')
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
subplot(2,1,2),[C,h] = contour(Var,10);
clabel(C,h);
xlabel('Ordrede la partie AR: na')
ylabel('Ordre de la partie exogène: nb')
legend('Variance de l erreur de prédiction')

% Var;  %// afficher la variance
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************

%*************************************************************************
%//////////////////////////////////////////////////////////////////////////
function [P,Q]=Charger

load('Peq_nonglacier_gletsch_horaire_90_95.txt');  
Peq_nonglacier_gletsch=Peq_nonglacier_gletsch_horaire_90_95;


%// Charger les Pluie equivalente de la partie glaciere (pluie +fonte neige)
load('Peq_glacier_gletsch_horaire_90_95.txt');  
Peq_glacier_gletsch=Peq_glacier_gletsch_horaire_90_95;




P=(Peq_nonglacier_gletsch+Peq_glacier_gletsch);

load('Qobs_gletsch_horaire_90_95.txt'); %Débits à l'exutoire des deux bassins (Gletsh/Alptal)
Q=Qobs_gletsch_horaire_90_95;  %Débits du Bassin (1)Alptal, (2)Gletsch


figure,subplot(2,1,1), plot(Q), 
ylabel('Debits à l exutoire du Bassin Gletsch'), grid on
subplot(2,1,2), plot(P)
ylabel('Precipitation Effective'), grid on
xlabel('Temps en pas d1 jour')
%/////////////////////////////////////////////////////////////////////////////
%**************************************************************************
%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
function [Teta]=ConditionInitiale(na,nb,nk,N1,Q,P);


for j=N1(1):N1(2)
    
    Y(j-N1(1)+1)=Q(j);

for k=1:na
Phi((j-N1(1)+1),k) =[-Q(j-k)];
end    
    

for k=1:nb
Phi((j-N1(1)+1),k+na) =[P(j-k+1-nk)];
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
%//////////////////////////////////////////////////////////////////////////
function [yh]=Simulation(N1,P,Q,na,nb,nk,Teta)
      

for i=(N1(1)+1):N1(2)   

%/////// Estimation de Phi //////////////////

for k=1:na        
Ph(k) =[-Q(i-k)]  ;  
end

for k=1:nb
Ph(k+na) =[P(i-k-nk+1)];
end
 

%/// Estimation des erreurs de prévisions du modele ARX
yh(i)=Ph*Teta;   %// Prévision

end

%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************


%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
function [VarE,AIC]=Critere(yh,Q,N1,na,nb);

%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%// Residus e= Q-yh   ///////////////////////////////////////////////////
%////////////////////////////////////////////////////////////////////////
e= [Q(N1(1)+1: N1(2))]- [yh(N1(1)+1:N1(2))'] ;
E=sum(e);
N=length(e);
E1=(e'*e);
%/// Variance des erreurs  ////////////////////
VarE=E1/N;
AIC=log10(E)+2*((na+nb)/N);
%figure, plot([N1(1)+1: N1(2)],Q(N1(1)+1: N1(2)));
%hold on, plot([N1(1)+1: N1(2)], [yh(N1(1)+1:N1(2))'], 'r')
%/////////////////////////////////////////////////////////////////////////////
%*****************************************************************************
