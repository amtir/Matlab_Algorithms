


function Degradation_Erreur_Pred

%**************************************************************************
%**************************************************************************
%/////  Initialisation /////
%// Choix d'une période de Crue 
N1=[2773,2950];  %// Crue du 26 avril 1994
%%// Coefficient de pondération
Lamda=1;

%// Choix de l'ordre du modèle ARMAX
na=2;
nb=2;
nc=0;
nk=1;
horiz=20;
%**************************************************************************
%**************************************************************************

%//  Programme principale

%// Prévision de 72h à partir de l'instant N1(1)
[e,yh,Q,u,U,St]=MMCRP([N1(1),N1(1)+horiz],na,nb,nk,Lamda);
yh=yh./10;
Q=Q./10;
N=length(e);
E(N,7)=0;
S(N,7)=0;
E(:,1)=e;
S(:,1)=u;
S2(:,1)=U./10;

figure(10), grid on, 
plot((N1(1)-10: N1(2)),Q(N1(1)-10: N1(2)));
hold on, plot(N1(1):N1(1)+horiz,yh(N1(1):N1(1)+horiz),'k');


%hold on, plot(N1(1):N1(1)+horiz,[yh(N1(1):N1(1)+horiz)+St(N1(1):N1(1)+horiz)],'c');
%hold on, plot(N1(1):N1(1)+horiz,[yh(N1(1):N1(1)+horiz)-St(N1(1):N1(1)+horiz)],'g');

title('Comparaison entre le débit mesuré et le débit prédit pour une prévision de 72h');
xlabel('Temps: Pas de temps horaire');
ylabel('Débit (mm) ');



%// Prévision de 72h à partir de l'instant N1(1)+i

for j=1:10
    
    i=j;
    
[e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+horiz+i],na,nb,nk,Lamda);
E(:,j)=e;
S(:,j)=u;
S2(:,j)=U./10;
yh=yh./10;

figure(10);
hold on, plot(N1(1)+i:N1(1)+i+horiz,yh(N1(1)+i:N1(1)+i+horiz),'k');

end
% 
% %// Prévision de 72h à partir de l'instant N1(1)+i
% i=3;
% 
% [e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+72+i],na,nb,nk,Lamda);
% yh=yh./10;
% Q=Q./10;
% E(:,2)=e;
% S(:,2)=u;
% S2(:,2)=U./10;
% 
% 
% figure(10);
% hold on, plot(N1(1)+i:N1(1)+i+72,yh(N1(1)+i:N1(1)+i+72),'g');
% 
% %// Prévision de 72h à partir de l'instant N1(1)+i
% i=6;
% 
% [e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+72+i],na,nb,nk,Lamda);
% yh=yh./10;
% Q=Q./10;
% E(:,3)=e;
% S(:,3)=u;
% S2(:,3)=U./10;
% 
% figure(10);
% hold on, plot(N1(1)+i:N1(1)+i+72,yh(N1(1)+i:N1(1)+i+72),'k');
% 
% %// Prévision de 72h à partir de l'instant N1(1)+i
% i=12;
% 
% [e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+72+i],na,nb,nk,Lamda);
% yh=yh./10;
% Q=Q./10;
% E(:,4)=e;
% S(:,4)=u;
% S2(:,4)=U./10;
% 
% figure(10);
% hold on, plot(N1(1)+i:N1(1)+i+72,yh(N1(1)+i:N1(1)+i+72),'y');
% 
% %// Prévision de 72h à partir de l'instant N1(1)+i
% i=18;
% [e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+72+i],na,nb,nk,Lamda);
% yh=yh./10;
% Q=Q./10;
% E(:,5)=e;
% S(:,5)=u;
% S2(:,5)=U./10;
% 
% figure(10);
% hold on, plot(N1(1)+i:N1(1)+i+72,yh(N1(1)+i:N1(1)+i+72),'c');
% 
% %// Prévision de 72h à partir de l'instant N1(1)+i
% i=24;
% [e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+72+i],na,nb,nk,Lamda);
% yh=yh./10;
% Q=Q./10;
% E(:,6)=e;
% S(:,6)=u;
% S2(:,6)=U./10;
% 
% figure(10);
% hold on, plot(N1(1)+i:N1(1)+i+72,yh(N1(1)+i:N1(1)+i+72),'k');
% 
% %// Prévision de 72h à partir de l'instant N1(1)+i
% i=36;
% [e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+72+i],na,nb,nk,Lamda);
% yh=yh./10;
% Q=Q./10;
% E(:,7)=e;
% S(:,7)=u;
% S2(:,7)=U./10;
% 
% figure(10);
% hold on, plot(N1(1)+i:N1(1)+i+72,yh(N1(1)+i:N1(1)+i+72),'k');
% 
% %// Prévision de 72h à partir de l'instant N1(1)+i
% %i=72;
% %[e,yh,Q,u,U,St]=MMCRP([N1(1)+i,N1(1)+72+i],na,nb,nk,Lamda);
% %yh=yh./10;
% %Q=Q./10;
% %E(:,8)=e;
% %S(:,8)=u;
% %S2(:,8)=U./10;
% 
% %figure(10);
% %hold on, plot(N1(1)+i:N1(1)+i+72,yh(N1(1)+i:N1(1)+i+72),'k');



%///  Plot des Erreurs    /////

figure;
title('Plot des erreurs en mm de prédiction pour des prévisions de 72h');
xlabel('Temps: Pas de temps horaire');
ylabel('Erreur de prédiction en mm');
for i=1:7 %8
  hold on,  plot(E(:,i));
end

%///  Plot des Erreurs    /////
figure;
title('Plot des erreurs relatives en % de prédiction pour des prévisions de 72h');
xlabel('Temps: Pas de temps horaire');
ylabel('Erreur relative de prédiction en %');
for i=1:7 %8
  hold on,  plot(S(:,i));
end


figure, boxplot(E);
title('Boxplot des erreurs de prédiction pour des prévisions de 72h ');
xlabel('Décalage entre le début de chaque prévision: 0, 3, 6, 12, 18, 24, 36, 72 heures');
ylabel('Boxplot pour chaque prévision');
figure, boxplot(E')
title('Boxplot des erreurs de prédiction pour différents horizon de prévision: 1--> 20h');
xlabel('Temps: Pas de temps horaire');
ylabel('Boxplot pour chaque pas horaire ');

figure, boxplot(S);
xlabel('Décalage entre le début de chaque prévision: 0, 3, 6, 12, 18, 24, 36, 72 heures');
figure, boxplot(S');


%///////// Boxplot du cumul des Erreurs relatives //////////////////////////////
figure, boxplot(S2);
xlabel('Décalage entre le début de chaque prévision: 0, 3, 6, 12, 18, 24, 36, 72 heures');
figure, boxplot(S2');
