function [F,P,V,U]=roughmovie(f1,f2,t);

%[f,x,y]=rsgeng2D(400,400,rmsr,20,20);
l=size(f1);
surf(f1+f2)
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
sc=5e9;
loops = t;
F(loops) = struct('cdata',[],'colormap',[]);
Z1=f1;
Z2=f2;
for y=1:l(2);
     X(:,y)=linspace(0,l(1)*5e-5,l(1));  %Sets up grid for x component
end

for x=1:l(1);
    Y(x,:)=linspace(0,l(2)*5e-5,l(2));  %Sets up grid for y component
end

for jj = 1:loops;
%     for k=1:l(2);
%         Z1new(:,k)=zeros(1,l(1));
%         Z2new(:,k)=Z1new(:,k);
%         Z1new(1,k)=Z1(end,k);
%         Z1new(2:end,k)=Z1(1:end-1,k);
%         Z2new(:,k)=Z2(:,k);
%         Z2new(end,k)=Z2(1,k);
%     end
%     Z1=Z1new;
%     Z2=Z2new;
     data=Z1+Z2;
     ZZ=data;
    ac=1-(30e6/5e9);
    %ac=1-(0.1+(jj-1)*0.005);
    ZZ=reshape(ZZ,l(1)*l(2),1);
    [FF,XX]=ecdf(ZZ);
    n=1;
    for ii=1:numel(FF);
        if FF(ii)>=ac; %Extract points in top percentage contact 
            cutt(n)=XX(ii); %assign values of height greater than asperity contact height
            n=n+1; %continue counter
        else
        end
    end
    cutt=min(cutt);
    data1=data;
    data(data<cutt)=0;
    data(data>=cutt)=1;
    data1(data1<cutt)=nan;
    du=data1-cutt;
    sigma=du.*60e9;
    %h = waitbar(0,'Grid progress');
rold=0;
    for i=1:l(1);
        for j=1:l(2);
           X1=X-(i*5e-5);   
           Y1=Y-(j*5e-5);
           dista=sqrt((X1.^2)+(Y1.^2)).*data; 
           %dir=(atan(Y1./X1));  %Add-in to establish direction to nearest neighbour
           dista(dista==0)=nan;
           r(i,j)=min(min(dista));
           rnew=r(i,j);
           if rnew>rold;
               rold=rnew;
               Xr=i;
               Yr=j;
           else
           end
        end
    %f=(i*j)/numel(data);
    %waitbar(f)
    end
    
   rr=max(max(r));
   V(jj)=rr;
   xr=linspace(0,l(1),10000);
   yr=((rr*(1./5e-5))^2-(xr-Yr).^2).^(1/2)-Xr;
   yr(real(yr)==-Xr)=nan;
   %close(h);
    figure(1)
    subplot(1,2,1)
    imagesc(sigma)
    hold on
    plot(xr,-real(yr),xr,real(yr)+2*Xr);
    shading flat; axis equal
    xlim([0 l(2)])
    ylim([0 l(1)])
    title(['\sigma_n=', num2str(((1-ac).*sc)/1e6) ,'MPa']);%, slip=',num2str(jj*50),'\mu m']);
    hold off
    colorbar
    
    %subplot(1,3,2)
%     pcolor(r)
   % scatter(Yr,Xr,200)
   % shading flat; axis equal
    %title(['\sigma_n=', num2str(((1-ac).*sc)/1e6) ,'MPa, slip=',num2str(jj*50),'\mu m']);
    %xlim([0 l(2)])
    %ylim([0 l(1)])
    figure(1)
    subplot(1,2,2)
    title(['Void scaling'])
    %scatter(jj*50,rr);
    scatter((((1-ac).*sc)/1e6),rr);
    %xlabel(['Displacement [\mu m]']);
    xlabel(['Normal stress [MPa]']);
    ylabel(['Void length [m]']);
    %ylim([10^(log10(V(1))-0.5) 10^(log10(V(1))+0.5)]);
    ylim([1e-4 1e-2])
    xlim([0 ((0.001+(loops-1)*0.005)*sc)/1e6])
    hold on
    drawnow
    F(jj) = getframe(gcf);
    P(jj)=(1-ac).*sc;
    U(jj)=jj*5e-5;
    figure(2)
    imagesc(data)
    hold on
    plot(xr,-real(yr),'r',xr,real(yr)+2*Xr,'r');
    shading flat; axis equal
    xlim([0 l(2)])
    ylim([0 l(1)])
    title(['\sigma_n=', num2str(((1-ac).*sc)/1e6) ,'MPa']);%, slip=',num2str(jj*50),'\mu m']);
    hold off
    colorbar
end

%     for j=1:l(1);
%         Z1new(j,:)=zeros(l(2),1);
%         Z2new(j,:)=Z1new(j,:);
%         Z1new(j,1:clip)=Z1(j,end-999:end);
%         Z1new(j,1001:end)=Z1(j,1:end-1000);
%         Z2new(j,1:end-1000)=Z2(j,1001:end);
%         Z2new(j,end-999:end)=Z2(j,1:1000);
%     end
%     Z1=Z1new;
%     Z2=Z2new;
%     data=Z1+Z2;
%     ZZ=data;
%     topheight=1-(ac); %Get percentage area contact of sample based on normal stress 
%     zz=reshape(data,l(1)*l(2),1);%Reshape data to perform CDF
%     [FF,XX]=ecdf(zz); %Run Kaplan-Maeir survivor CDF estimate
%     n=1; %Initialise counter for asperity extraction 
%     for ii=1:numel(FF);
%         if FF(ii)>=topheight; %Extract points in top percentage contact 
%             cutt(n)=XX(ii); %assign values of height greater than asperity contact height
%             n=n+1; %continue counter
%         else
%         end
%     end
%     cutt=min(cutt); %Extract minimum of apserity heights in distribution
%     data(data>cutt) = 2; %Assign 2's to asperity start
%     data(data<=cutt) = 1;