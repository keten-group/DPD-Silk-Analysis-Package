close all;
load('boundary_condition.txt');
Lx=boundary_condition(1);Lx0=Lx;
Ly=boundary_condition(2);Ly0=Ly;
Lz=boundary_condition(3);Lz0=Lz;
L=[Lx Ly Lz];

load('cluster_connectivity.txt');
load('cluster_coordinate.txt');
ADJ=cluster_connectivity>0;[imax,jmax]=size(ADJ);
XY=cluster_coordinate;
additional_XY_orig=[];
additional_XY_pbc =[];
count=0;count1=0;multi=zeros(1,1000);

for i=1:imax
    for j=i+1:jmax
        dr=XY(i,:)-XY(j,:);
        if ((ADJ(i,j)>0) && max((max(abs(dr))> 0.5*L))>0 )
            count=count+2;
            multi(count-1)=cluster_connectivity(i,j);
            multi(count)=cluster_connectivity(i,j);
            ADJ(i,j)=0;
            ADJ(j,i)=0;           
            additional_XY_orig=[additional_XY_orig;XY(i,:)];
            additional_XY_orig=[additional_XY_orig;XY(j,:)];
            sign_i=sign(XY(i,:));
            sign_j=sign(XY(j,:));
            %sign_ij=(sign_i.*sign_j)<0;
            sign_ij=(abs(dr))> 0.5*L;
            XYj_p=XY(j,:)-sign_ij.*sign_j.*L;
            XYi_p=XY(i,:)-sign_ij.*sign_i.*L;
            additional_XY_pbc =[additional_XY_pbc ;XYj_p];
            additional_XY_pbc =[additional_XY_pbc ;XYi_p];
        end
    end
end
multi=multi(1:count);
addX=zeros(2,count,3);
addX(1,:,:)=additional_XY_orig;
addX(2,:,:)=additional_XY_pbc;
[i, j] = find(triu(ADJ)); 
X = permute(cat(3, XY(i, :),XY(j, :)), [3 1 2]);

figure1=figure(1);
set(figure1,'Position',[100 500 1800 400]);
for index=1:nnz(ADJ)/2
ii=i(index);jj=j(index);
%cluster_connectivity(ii,jj)
plot3(X(:,  index, 1), X(:, index, 2), X(:, index, 3),'k-',...
    'LineWidth',0.5+0.25*(cluster_connectivity(ii,jj)-1));
hold on
end
for index=1:count
plot3(addX(:,  index, 1), addX(:, index, 2), addX(:, index, 3),'k-',...
    'LineWidth',0.5+0.25*(multi(index)-1));
hold on
end

plot3(XY(:,1),XY(:,2),XY(:,3),'ko','markerfacecolor','r')
hold on
boundary=[ -Lx/2 -Ly/2 0 ;
           -Lx/2  Ly/2*0.999 0 ;
            Lx/2  Ly/2*0.999 0 ;
            Lx/2 -Ly/2 0 ;
           -Lx/2 -Ly/2 0 ];
%plot3(boundary(:,1),boundary(:,2),boundary(:,3),'k:','LineWidth',0.5);
set(gca,'XTick',[],'YTick',[],'ZTick',[],...
    'XColor',[1 1 1],'YColor',[1 1 1 ],'ZColor',[1 1 1],...
    'Position',[0.5-Lx/Lx0/3/2*0.95 0.5-Ly/Ly0/2*0.95 Lx/Lx0/3*0.95 Ly/Ly0*0.95]);
%set(gca,'XTick',[],'YTick',[],'ZTick',[],...
%    'XColor',[1 1 1],'YColor',[1 1 1 ],'ZColor',[1 1 1]);
set(gcf,'PaperPosition',[0 0 27 6]);
hold off
xlim([-Lx/2 Lx/2*1.004]);
ylim([-Ly/2 Ly/2*1.004]);
zlim([-Lz/2 Lz/2]);
view(2);
daspect([1 1 1]);
print(figure1,'-dpng','movie.png');
close(figure1);
