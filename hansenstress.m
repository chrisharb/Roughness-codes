function sol = hansenstress(data,dx,SN,SN1,SN2);
%% Hansen stress script
%This code solves the elastic contact problem in 2d using the fast fourier
%transform algorith suggested by Stanley & Kato, and implemented by Hansen
%et al. 2002. The code calculates stress interference using a Green's
%function transformed into fourier space, multiplied by a unit
%displacement at points corresponding to asperities. The force required to
%closure acts as the multiplier for further surface deformations to give a
%new surface input, which is then resolved to find the new position of an
%asperity. Some features are still under development, and it is hoped to
%extend the model to work on real fault topography data sets.
%
%Inputs:
%       data = a matrix of evenly spaced elevation points in metres
%       dx = the sample spacing of the topography data
%       SN = the normal stress that the code should run to in Pa, do not 
%       use MPa otherwise the run will be rather short
%       SN1 = first normal stress for box counting, in Pa
%       SN2 = second normal stress for box counting, in Pa
%Outputs:
%       r = the maximum distance between asperities which is equivalent to 
%       lambda c in Harbord et al. (2017)
%       minSc = the stress multiplier used for each iteration
%       Sc{} = the minimum multiplier for contact at each point at every
%       iteration
%       P = the stress distribution at each iteration
%       Znew = the topography at each iteration, which can be used to make
%       a 3d visulisation of surface topography with increasing load
%       vx = the x coordinates of voronoi vertices at each iteration
%       vy = the y coordinates of voronoi vertices at each iteration
%       Pdf = a probablity density function of asperities- may be used in
%       conjunction with Beeler et al. 2008 model of Flash weakening to
%       look at critical weakening velocities
%       Pdfe = bin edges of Pdf
%       Sn = the average normal stress across the entire surface, which
%       gives gross normal stress
%       closure = the normal component of surface closure, to check
%       results against those of Yoshioka & Scholz, also in development...
%       
%Optional outputs dependant whether the box counting algorith has been
%triggered by the user:
%       pressC = probability density lines to give density of data
%       pressa = average stress within each box count sub-domain
%       minp = the minimum stress within each box count sub-domain
%       maxp = the maximum stress within each box count sub-domain
%       cella = the linear dimension of each sub-domain, for plotting of
%       homogenisation graphs
%       suba = the number of asperities contained within each subdomain 
%       length, for identification of a homogenisation scale
%% Begin the code
tic
Z = data;
Zmin = min(min(Z));
Zmax = max(max(Z));
Z = Z-Zmin;
% figure(2)
% surf(Z)
% shading flat
g = 50e9; %Youngs modulus of substrate
cut=3; %Number of points to take from top
l = size(Z); %Get the size of input data
%Reshape the data to get a 1d list for asperity picking
Zlist = reshape(Z',l(1)*l(2),1);
%Sort list in ascending order, storing positions of each point
[Zsort,ZsortI] = sort(Zlist); 
Zmax = Zsort(end-cut+1:end); %Extract values of elevation within cutoff
%Extract corresponding positions through 1-D array
ZmaxI = ZsortI(end-cut+1:end);
%Calculate j indexes of points by remainder of index division
Zmaxi = rem(ZmaxI,l(1));
Zmaxi(Zmaxi==0)=l(1); %If no remainder j = l(1) as multiple of l(1)
Zmaxj = ((ZmaxI-Zmaxi)./l(1))+1; %Calculate i indexes

Zmaxi=Zmaxi';
Zmaxj=Zmaxj';
Area = l(1)*dx*l(2);
SN30 = 0;
SN100 = 0;

%% (Un)comment this section to (not) initialise a movie
% ax = gca;
% ax.NextPlot = 'replaceChildren';
%F(loops-1) = struct('cdata',[],'colormap',[]);


%% Stress kernels a la Hansen et al. 2002 for surface stress and 
% displacement calculations

for i = 1:l(1);
    for j = 1:l(2); %Caluculate greens function kernel through mesh
        G(i,j) = sqrt(((i-1)*dx)^2+((j-1)*dx)^2)^-1+...
            sqrt(((i-1)*dx)^2+((l(2)+1-j)*dx)^2)^-1+...
            sqrt(((l(1)+1-i)*dx)^2+((j-1)*dx)^2)^-1+...
            sqrt(((l(1)+1-i)*dx)^2+((l(2)+1-j)*dx)^2)^-1; 
    end
end

G(1,1)=4; %To remove singularity in grid
G=((1-0.25^2)/pi*g).*G; %Stress correction factor
FG = fft2(G); %2-d FFT of kernel for stress calculations later on

%% Displacement matrix

U = zeros(size(G)); %Create solution array for displacement matrix
du = 1e-8; %Assign unit displacement to initial asperity location
Znew{1} = Z; %Initialise elevation cell matrices
P{1} = zeros(size(Z)); %Initialise pressure cell-array

for y=1:l(2);
     X(:,y)=linspace(0,l(1)*dx,l(1));  %Sets up grid for x component
end

for x=1:l(1);
    Y(x,:)=linspace(0,l(2)*dx,l(2));  %Sets up grid for y component
end
%reshape grid of X coordinates into 1-d array for calculations later
Xr = reshape(X,l(1)*l(2),1);  
%reshape grid of Y coordinates into 1-d array for later calculations
Yr = reshape(Y,l(1)*l(2),1);   

h=waitbar(0,'[Insert witty title here]');   %Initialise waitbar

f=1; %Initialise counter
normalstress = 0; %Initialise normal stress tracker
closure = 0; %Initialise cloasure tracker
%% Run the loading loop 

while normalstress < SN;
    if normalstress>SN1; %To run box counter at normal stress of interest
        SN30 = SN30+1;
    else
    end
    
    if normalstress>SN2;%To run box counter at normal stress of interest
        SN100 = SN100+1;
    else
    end
    
    f=f+1; %advance the counter
    %Update the wait bar
    waitbar(normalstress/SN,h,['At \sigma_n = ' ...
        num2str(normalstress/1e6,3) ' MPa of ' num2str(SN/1e6,3)...
        ' MPa, time elapsed = ' num2str(toc,3)])
    Uu{f} = U;

    for p = 1:length(Zmaxi);
        ii = Zmaxi(p);  %Extract i indexes of asperities
        jj = Zmaxj(p);  %Extract j indexes of asperities
        Uu{f}(ii,jj)=du; %Impose displacement change at each asperity
    end
   
    FUU = fft2(Uu{f}); %2d FFT of displacement matrix
    Fhe = FG.*FUU;  %Multiply kernel FFT by displacement matrix FFT
    P{f} = real(ifft2(Fhe));  %2d inverse FFT of force matrix FFT
    P3 = P{f}; %Copy pressure map to new array for other operations
    %To exclude asperity stresses from minimum point calculation as will 
    %grow singularity
    for p = 1:length(Zmaxi); 
        ii = Zmaxi(p);  %Extract i indexes of asperities
        jj = Zmaxj(p);  %Extract j indexes of asperities
        P{f}(ii,jj)=nan; %Replace asperity stress with NaN to ignore
        P3(ii,jj)=5e9; 
    end
    %Uncomment below If you want to visualise stress distribution at every 
    %step
    %surf(P{f}); shading flat; daspect([1 1 1e7]); 
    %zlim([0 1e8]);xlim([0 l(2)]);ylim([0 l(1)]); colorbar;
    %% Minimum stress solution
    %Minimum multiplier at each point to form new asperity
    Sc{f} = ((P{f})./(g)).*(Zmax(1)-Znew{f-1}); 
    %Extract minimum value of stress multiplier for the new asperity
    minSc(f) = min(min(Sc{f})); 
    %Calculate the closure (this needs checking in future)
    closure(f+1) = minSc(f).*du;
    %Find the matrix position of the new asperity
    [ScI,ScJ] = find(Sc{f}==minSc(f)); 
    %Calculate the new surface topography given the applied stress state
    Znew{f} = Znew{f-1}+((minSc(f).*(P{f}./g)));
    Zmaxi = horzcat(Zmaxi,ScI);
    Zmaxi2 = Zmaxi.*dx;
    Zmaxj = horzcat(Zmaxj,ScJ);
    Zmaxj2 = Zmaxj.*dx;
    

%% Geometric calculations using MATLAB inbuilt voronoi function
clear vx vy
%run the voronoi function to extract the maximum distance locations
[vx,vy] = voronoi(Zmaxi.*dx, Zmaxj.*dx); 
fff = size(vx); %extract size of the voronoi point array
%reshape the voronoi vertices for further processing
vx=reshape(vx,fff(1)*fff(2),1); 
%reshape the voronoi vertices for further processing
vy=reshape(vy,fff(1)*fff(2),1); 
%Loop to remove voronoi vertices outside of the range of the data
for i = 1:fff(1)*fff(2);
    if vx(i) > (l(1)*dx) || vx(i)<dx || vy(i)<dx || vy(i)>(l(2)*dx)
    vy2(i)=0;
    vx2(i)=0;
    else
    vx2(i)=vx(i);
    vy2(i)=vy(i);
    end
end
%Replace with empty after the loop
vy2(vy2==0)=[];
vx2(vx2==0)=[];

clear r dista
dista2=[];
dista3=[];
%Calculate distances from every asperity to voronoi vertices
for i = 1:numel(vy2);
    dista = sqrt((Zmaxi2-vx2(i)).^2+(Zmaxj2-vy2(i)).^2);
    r(i) = min(min(dista));
end

%Calculate distances for SN1
if SN30 == 1;
    for i = 1:numel(Zmaxi);
    dist = sqrt((Zmaxi2-Zmaxi2(i)).^2+(Zmaxj2-Zmaxj2(i)).^2);
    dista30(:,i) = dist;
    end
else
end

%Calculate distances for SN2
if SN100 == 1;
    for i = 1:numel(Zmaxi);
    dist = sqrt((Zmaxi2-Zmaxi2(i)).^2+(Zmaxj2-Zmaxj2(i)).^2);
    dista100(:,i) = dist;
    end
else
end
    
%% Box counting routines to establish a homogenisation scale    
    
% Pp = P{f};
% Pp(isnan(Pp))=1e9;
P2 = P{f};
P2(P2==5e9)=1;
P2(isnan(P2))=1;
sn(f) = mean(mean(P2));
normalstress = sn(f);
P2(P2>1)=0;
P2(P2<1)=0;

% Asperity density calculator
% Half spaces
P4 = P{f};
cell = 0.005; %Define initial cell dimensions
clear suba subci subcj

if SN30 == 1; %Run box counter 1 if required
    minp30 = [];
    maxp30 = [];
    suba30 = [];
    cella = [];
    pressa30 = [];
for i = 1:100;
    clear subasps ci cj
    X1 = round(cell/dx); %Round cell to nearest point in x dir
    Y1 = round(cell/dx); %Round cell to nearest point in y dir
    posI = [1:X1:l(1)];  %Create extraction vector for x dir
    posJ = [1:Y1:l(2)];  %Create extraction vector for y dir
    for iii = 2:numel(posI); %Run box calculator for each sub domain
        for jjj = 2:numel(posJ);
            %Calculate number of asperities in each sub domain by summation
            subasps(iii-1,jjj-1) =  sum(sum(P2(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj))))./cell.^2; 
            %Calculate average stress in each sub domain
            press(iii-1,jjj-1) =  mean(mean(P3(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj)))); 
            %Calculate maximum stress in each sub domain
            maxp(iii-1,jjj-1) = max(max(P3(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj)))); 
            %Calculate minimum stress in each sub domain
            minp(iii-1,jjj-1) = min(min(P3(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj)))); 
            %Define x positions at iii, jjj
            ci(iii-1,jjj-1) = posI(iii-1)+(posI(iii)-posI(iii-1))/2; 
            %Define y positions at iii, jjj
            cj(iii-1,jjj-1) = posJ(jjj-1)+(posJ(jjj)-posJ(jjj-1))/2; 
        end
    end
    %Probability intervals for 
    pressC(i,:) = prctile(reshape(press,numel(press),1),...
        [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],1);
     %suba{i} = reshape(subasps,numel(subasps),1);
     minp30 = horzcat(minp30,reshape(minp,1,numel(minp)));
     maxp30 = horzcat(maxp30,reshape(maxp,1,numel(maxp)));
     suba30 = horzcat(suba30,reshape(subasps,1,numel(subasps)));
     %suba30(suba30<1)=nan;
     subci30{i} =  reshape(ci,numel(ci),1);
     subcj30{i} =  reshape(cj,numel(cj),1);
     pressa30 =  horzcat(pressa30,reshape(press,1,numel(press)));
     cell2 = ones(numel(minp),1).*cell;
     cella = horzcat(cella,cell2'); 
%      figure(20+i);
%      scatter(subci{i},subcj{i},suba{i}*100,'o');
    cell = cell-5e-5; %Reduce cell size for next iteration
end
end

if SN100==1; %Run box counter 2 if required
    clear maxp minp press subasps ci cj
    minp100 = [];
    maxp100 = [];
    suba100 = [];
    pressa100 = [];
    cella2 = [];
for i = 1:100;
    clear subasps ci cj
    X1 = round(cell/dx); %Round cell to nearest point in x dir
    Y1 = round(cell/dx); %Round cell to nearest point in y dir
    posI = [1:X1:l(1)];  %Create extraction vector for x dir
    posJ = [1:Y1:l(2)];  %Create extraction vector for y dir
    for iii = 2:numel(posI); %Run box calculator for each sub domain
        for jjj = 2:numel(posJ);
            %Calculate number of asperities in each sub domain by summation
            subasps(iii-1,jjj-1) =  sum(sum(P2(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj))))./cell^2; 
            %Calculate average stress in each sub domain
            press(iii-1,jjj-1) =  mean(mean(P3(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj)))); 
            %Calculate maximum stress in each sub domanin
            maxp(iii-1,jjj-1) = max(max(P3(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj))));
            %Calculate minimum stress in each sub domain
            minp(iii-1,jjj-1) = min(min(P3(posI(iii-1):posI(iii),...
                posJ(jjj-1):posJ(jjj))));
            %Calculate x positions at iii, jjj
            ci(iii-1,jjj-1) = posI(iii-1)+(posI(iii)-posI(iii-1))/2;
            %Calculate y positions at iii, jjj
            cj(iii-1,jjj-1) = posJ(jjj-1)+(posJ(jjj)-posJ(jjj-1))/2; 
        end
    end
    pressC100(i,:) = prctile(reshape(press,numel(press),1),...
    [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],1);
     %suba{i} = reshape(subasps,numel(subasps),1);
     minp100 = horzcat(minp100,reshape(minp,1,numel(minp)));
     maxp100 = horzcat(maxp100,reshape(maxp,1,numel(maxp)));
     suba100 = horzcat(suba100,reshape(subasps,1,numel(subasps)));
     %suba30(suba30<1)=nan;
     subci100{i} =  reshape(ci,numel(ci),1);
     subcj100{i} =  reshape(cj,numel(cj),1);
     pressa100 =  horzcat(pressa100,reshape(press,1,numel(press)));
     cell2 = ones(numel(minp),1).*cell;
     cella2 = horzcat(cella2,cell2'); 
%      figure(20+i);
%      scatter(subci{i},subcj{i},suba{i}*100,'o');
    cell = cell-5e-5; %Reduce cell size for next iteration
    
end
end


rr(f) = max(r);
%% Uncomment this section to build a movie
%     hh = figure(11);
%     set(hh, 'Position', [100 300 666 500])
%     subplot(1,2,1)
%     imagesc(P2)
%     shading flat; axis equal
%     xlabel(['X Distance x' num2str(dx) ' m']);
%     ylabel(['Y Distance x' num2str(dx) ' m']);
%     xlim([0 l(2)])
%     ylim([0 l(1)])
%     title(['Surface with ' num2str(f+1) ' asperities, \sigma_n = ' num2str(sn(f)/1e6) ' MPa']) %, \Delta\sigma = ' num2str((max(max(P3))-min(min(P3)))/1e6) ' MPa' ]);%, slip=',num2str(jj*50),'\mu m']);
%     hold off
%     c=colorbar;
%     ylabel(c,'Stress [MPa]')
%     caxis([0 max(max(real(P{f})))])
% 
% %     subplot(1,3,2)
% %     voronoi(Zmaxj,Zmaxi)
% %     axis equal
% %     xlim([0 l(2)])
% %     ylim([0 l(1)])
% %     title(['Voronoi map with ' num2str(f) ' asperities']);
%     
%     subplot(1,2,2)
%     plot(sn./1e6, 2.*rr)
%     %ylim([0 1.5e-2])
%     %xlim([0 200])
%     daspect([200 2/3*1.5e-2 1])
%     ylabel('Void length scale [m]')
%     xlabel('Averaged normal stress [MPa]')
%     hold off
%     drawnow
%     F(f-1) = getframe(gcf);
%     
%%
PP = reshape(P{f},l(1)*l(2),1)/max(max(P{f}));
[N,edges] = histcounts(PP,'Normalization','probability');
Pdf{f} = N;
Pdfe{f} = edges;
end

close(h)
%% Plotting of results
%Plot asperity map
figure(1)
scatter(Zmaxi,Zmaxj,'o'); axis equal
xlim([0 l(1)]);
ylim([0 l(2)]);

%Plot final voronoi map
figure(2)
voronoi(Zmaxi,Zmaxj)

%Plot void length as function of the number of contacts
figure(3)
plot(rr)
ylabel('Void length')
xlabel('Number of contacts')

%Mesh the final stress state
figure(4)
surf(real(P{f}));
shading flat
zlabel('Normal stress [Pa]')
daspect([1 1 1e8])
xlim([0 l(2)]);
ylim([0 l(1)]);

% figure(6)
% surf(pressa{end})
% axis equal
% shading flat
% view([0 0 1])
% figure(5)
% plot(sn./1e6,cumsum(closure).*1e6)
% xlabel('Normal stress [MPa]')
% ylabel('Closure [\mum]')
% title('Closure of surface')

% dista30a = reshape(dista30,numel(dista30),1);
% dista100a = reshape(dista100,numel(dista100),1);
% [xx30,yy30] = ecdf(dista30a);
% [xx100,yy100] = ecdf(dista100a);
% figure(100)
% plot(yy30,xx30)
% hold on
% plot(yy100,xx100)
% ylabel('Cumulative probability')
% xlabel('Length scale (m)')
% legend('30 MPa','100 MPa')
%% Save results to solution array

sol.r=rr;
sol.minSc=minSc;
sol.Sc=Sc;
sol.P=P;
sol.Znew=Znew;
sol.vx=vx;
sol.vy=vy;
sol.Pdf = Pdf;
sol.Pdfe = Pdfe;
%ans.F = F;
sol.sn=sn;
sol.closure=cumsum(closure);
%sol.dista30 = dista30;
%sol.dista100 = dista100;

if SN > SN30;
sol.pressC = pressC;
sol.pressa30 = pressa30;
sol.minp30 = minp30;
sol.maxp30 = maxp30;
sol.cella = cella;
sol.suba30 = suba30;
elseif SN > SN100;
sol.pressa100 = pressa100;
sol.minp100 = minp100;
sol.maxp100 = maxp100;
sol.pressC100 = pressC100;
sol.cella2 = cella2;
sol.suba100 = suba100;
figure(101)
title('30 & 100 MPa')
% scatter(cella,minp30)
% hold on
% scatter(cella,maxp30)
scatter(cella,pressa30)
hold on
%scatter(cella2,pressa100)
else
end
toc

end  
