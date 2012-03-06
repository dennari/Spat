
offsetX = 3478085;
offsetY = 7085263;
bear2010 = csvread('BearObs2010XY.csv');
bear2010 = bear2010-repmat([offsetX offsetY],size(bear2010,1),1);
%plot(bear2010(:,1)+offsetX,bear2010(:,2)+offsetY,'.','MarkerSize',3)
%hold on;
%draw_prov_borders;
%error('sdfsdf');
gridn = 20;
[P,PQ,XT] = lgcp(bear2010,'latent_method','Laplace','gridn',gridn);
w = sqrt(size(P,1));
xt1 = reshape(XT(:,1),w,w)+offsetX;
xt2 = reshape(XT(:,2),w,w)+offsetY;
P = reshape(P,w,w);
pcolor(xt1,xt2,P);
shading flat
colormap('jet')
cx=caxis;
cx(1)=0;
caxis(0.7^5*cx);
%colorbar
hold on;
draw_prov_borders;
%hold on;
%plot(bear2010(:,1)+offsetX,bear2010(:,2)+offsetY,'.r','MarkerSize',3)