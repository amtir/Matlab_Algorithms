
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%K-Plus_Proches_Voisins:  M�thode de pr�vision Stochastique non param�trique
%//////////////////////////////////////////////////////////////////////////
%Arguments d'entr�es:
%/////////  M1,M2,M3 sont l'ordre du syst�me //////////////////////////////
%////////   L repr�sente le temps de r�ponse de la neige fondue //////////
% ///////   K est le nombre des voisins    //////////////////////////////
% ///////   n represente le nombre de jour dans l'historique  //////////
% ////////  m le nombre de pr�vision journali�re  ///////////////////////
%//////////////////////////////////////////////////////////////////////////
%En sortie le vecteur pr�vision: Prev
%// Affichage sous forme graphique
%// Des crit�res de validation s'affichent dans la commande Window
%////////////////////////////////////////////////////////////////////
%Exemple: Pour     K_Plus_Proches_Voisins (M1,M2,M3,L,K,n,m)
%                  K_Plus_Proches_Voisins(2,2,0,0,3,3200,100)
%*************************************************************************
%*************************************************************************
function K_Plus_Proches_Voisins (M1,M2,M3,L,K,n,m)

[P,Q]=Charger(n,m);
H = waitbar(0,'Please wait...');
[Prev,h]=Algorithme (n,m,Q,P,K,M1,M2,M3,L);
Dessiner(P,Q,Prev,n,m);
Critere_Validation(Q,Prev,n,m,M1,M2,M3);
close(H) 

%**************************************************************************
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%///// Charger les donn�es  ///////////////////
%//////////////////////////////////////////////////////////////////////////
%Commencer par charger les donn�es: moyenne journali�re du bassin
function [P,Q]=Charger(n,m)


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
ylabel('Precipitation Equivalente (mm)'), grid on
xlabel('Temps en pas d1 jour: De Janvier 1992 au 31 d�cembre 2001')

%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%**************************************************************************
%**************************************************************************

function [W,h]=Algorithme (n,m,Q,P,K,M1,M2,M3,L)
g=0; %// Cette variable est utilise pour montrer la 

%Debut de notre Boucle pour les m previsions futures
for h= n+1:n+m   %// n+1 est la premiere prevision et n+m la derniere
g=g+1;  %// permet d'afficher la barre de progression
[V]=Construction_Vecteur(Q,P,M1,M2,M3,L,n,m,h); 
[Nombr_An,E,fen]=Nombre_Annee(h);
[E]=Recherch_Historiq(V,E,Nombr_An,fen,h,Q,P,M1,M2,M3);
[w]=Moyenne_K_Prem_Voisins(E,Q,K,g,m,h);
W(h)=w;
end

%///////////////////////////////////////////////////////////////////////////////////
%**************************************************************************
%**************************************************************************
%///////////////////////////////////////////////////////////////////////////////////
%Construction d'un vecteur caract�risant le comportement du Bassin 
%//////   versant � l'instant n 
%///////////////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
function [V]=Construction_Vecteur(Q,P,M1,M2,M3,L,n,m,h)
for i= 1:M1   
D(i)=Q(h-i);  %Vecteur qui contient les M1 d�bits pr�c�dent l'instant n+1
end

for i= 1:M2
Pr(i)= P(h-i);   %Vecteur qui contient les M2 pr�cipitations
end                 % Ici P(n+1) est une pr�vision meterologique 
 
for i = 1:M3
Te(i)=T(h-i-L); %Vecteur qui contient les M3 temperature avec un decalage de L
end        % Ici dans le cas ou L=0, Te(n+1) est une pr�vision meterologique      

%V=[D,Pr,Te];  % Vecteur qui r�sume l'�tat actuel (n) du Bassin
V=[D,Pr];
%V=[D];

%//////////////////////////////////////////////////////////////////////////
%**************************************************************************
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%//// Nombres d'annees dans le passee qui peuvent etre utilisees dans
%////  l'historique  ///////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
function [Nombr_An,E,fen]=Nombre_Annee(la_date)
fen=30;  %// fenetre de 30 jours (1 mois) � gauche et � droite
%// Nombre d'ann�e dans l'historique
for k=1:1000
 debut= la_date-(k*365)-fen;
 if ((debut <=10)&((la_date-((k-1)*365)-fen)>=10)|(debut==0))  %S'assurer que les donn�es existent
        Nombr_An=(k-1); %// Nombre d'ann�e dans l'historique
    end
 end  

 E(1:la_date)=100000000;  %// Construction et initialisation de E par d�faut un grand nombre

%**************************************************************************
%**************************************************************************
 %/////////////////////////////////////////////////////////////////////////
 %/////////  Recherche dans l'historique  /////////////////////////////////
 %/////////////////////////////////////////////////////////////////////////
 function [E]=Recherch_Historiq(X,E,Nombr_An,fen,la_date,Q,P,M1,M2,M3)
 
 for k=1:Nombr_An
  debut= la_date-k*365-fen;
  fin= la_date-k*365+fen;
 
for i =debut:fin  %Parcourt dans toute l'historique
    
    for j= 1:M1   
    D1(j)=Q(i-j);  %Vecteur qui contient les M1 d�bits pr�c�dent l'instant n
    end

    for j= 1:M2
    Pr1(j)= P(i-j);   %Vecteur qui contient les M2 pr�cipitations
    end
 
    for j = 1:M3
    Te1(j)=T(i-j-L); %Vecteur qui contient les M3 temperature avec un decalage de L
    end
        
        %H=[D1,Pr1,Te1];  % Vecteur qui r�sume l'�tat du Bassin � l'instant i
       H=[D1,Pr1];
        %H=[D1];
        
        d=X-H;  %Diff�rence entre le vecteur actuel(instant n) et le vecteur (instant i)
        d=d*d';  %Somme des diff�rences au 2 
        E(i)=sqrt(d) ;  %Calcul de la norme Euclidienne 
        
end

end

%**************************************************************************
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%////////recherche des K_voisins les plus proches//////////////////////////
%////////Comparaison des normes et mise dans l'ordre croissant des normes//
%/////// I :  Vecteur contenant l'ordre de croissance/////////////////////
%//////////////////////////////////////////////////////////////////////////
%////////////////////////////////////////////////////////////////////////
function [w]=Moyenne_K_Prem_Voisins(E,Q,K,g,m,h)

[C,I] = sort(E);   %// Ordonner le vecteur E, I contient l'ordre 
                   %// Voir le help dans matlab, la fonction sort

%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////Prevision du d�bit///////////////////////////
%///////////////// Moyenne des K_plus proches voisins ////////////////////////
%//////////////////////////////////////////////////////////////////////////

M=0;    % Initialisation � z�ro de la valeur pr�dite

%Calcul de la moyenne des K voisins les plus proches
for i=1:K     %Chercher les K voisins
    M=M+Q(I(i)+1);   %Sommer les K voisins : prendre Q[n+1]=Q[I(i)+ 1]
end                         % + 1 pour sommer le d�bit suivant!!! au voisin
M=M/K;   % Moyenne

%Stocker les Pr�visions dans la variable w
w=M;

waitbar(g/(m))     %/// permet d'afficher la barre de progression 

%//////////////////////////////////////////////////////////////////////////
%**************************************************************************
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%///  R�sultat sous forme graphique //////////////////////////////////////
%Tracer la courbe de pr�vision (1 jour) obtenue par la m�thode des K voisins les
%plus proches
function Dessiner(P,Q,w,n,m)

%/////// Plot des Pr�cipitations  //////////////
%figure,subplot(4,1,1),plot((n+1:n+m),P(n+1:n+m)); 
%title('Precipitation Equivalente ');
%ylabel('Precipitation Equivalente'), grid on
%xlabel('Temps en pas d1 jour')
%legend('Precipitation Equivalente (mm)') ;

figure,subplot(3,1,1), plot((n+1:n+m),Q(n+1:n+m)), grid on;
hold on, plot((n+1:n+m),w(n+1:n+m),'r')   %'--r');
%hold on, plot((n+1:n+m),w(n+1:n+m),'ko'); 
%title('Comparaison entre le d�bit mesur� et le d�bit pr�dit');
legend('D�bit mesur�','D�bit pr�dit') ;
xlabel('Temps: Pas de temps journalier');
ylabel('D�bit (mm)');

%/////// Calcul et plot des erreurs  //////////////
e(n+1:n+m)=Q(n+1:n+m)-(w(n+1:n+m))';
subplot(3,1,2),plot((n+1:n+m),100*(e(n+1:n+m))/max(Q(1:n+m))); 
%title('Erreur de pr�diction en %');
ylabel('Erreur de pr�diction en %'), grid on
xlabel('Temps en pas d1 jour')
legend('Erreur de pr�diction en %') ;
%//////////////////////////////////////////////////////////////////////////
%// Plot du debit mesure en fonction du debit simule
%//////////////////////////////////////////////////////////////////////////
subplot(3,1,3), plot(Q(n+1:n+m),w(n+1:n+m),'o');
f=max(Q);
%axis([0 f 0 f]);
%title('Relation entre le d�bit mesur� et le d�bit pr�dit');
legend('nuage de points') ;
xlabel('D�bit mesur� [mm]');
ylabel('D�bit pr�dit [mm]');


%**************************************************************************
%**************************************************************************
%//////////////////////////////////////////////////////////////////////////
%/////////  Calcul et affichache des criteres de validation dans la /////
%////////   commande window de matlab  //////////////////////////////////
function Critere_Validation(Q,w,n,m,M1,M2,M3)   %///////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////

%/////// Calcul et plot des erreurs  //////////////
e(n+1:n+m)=Q(n+1:n+m)-(w(n+1:n+m))';
%//////// calcul de la variance de l'erreur /////////////////////////////
N=length(e(n+1:end));
E1=(e(n+1:end)*(e(n+1:end))');
VarE=E1/N

%/// Le biais  //////////////////////////////////////////////////
m1=mean(Q(n+1:n+m));
m2=mean(w(n+1:n+m));
Biais=100*[(m2-m1)/m1]

%/////////////////////////////////////////////////////////////////////////
%// Coefficient de correlation 
%/////////////////////////////////////////////////////////////////////////
r=corrcoef(Q(n+1:n+m),w(n+1:n+m));
r=r(1,2)
%//////////////////////////////////////////////////////////////////////////

%///////////////////////////////////////////////////////////
%/// Coefficient de Nash-Suttcliffe  ///////////////////////
 E1=(e*e');
    E2=(Q(n+1:n+m));
E2=E2'*E2;
Nash= 1-(E1/E2) 

%/// AIC: Akaike Information Criterion  ////////////////////
E=sum(e);
E=abs(E);
AIC=log10(E)+2*((M1+M2+M3)/N)

%/// FPE: Final Prediction Error  ////////////////////
FPE=E1*((N+M1+M2+M3-1)/(N-M1-M2-M3-1))
%//////////////////////////////////////////////////////////////////////////
%*****************************************************************************
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////