function ans=maxellipse6(data1,dx,sc);

%Inputs:
%       data = a matrix comprised of areas of slip patch indicated by a 1,
%              and non-slipped patches by 0.
%       ar   = the aspect ratio of the ellipse used to define the maximum 
%              area between slipped patches. This is fundamentally
%              controlled by fracture mechanics. This feature will be added
%              in due course
%       dx   = the spacing of each data point in m
%       perarea = Area of real contact as a percentage of the nominal
%               contact area
%
%This code is being developed to establish the maximum spacing between
%points of slip patches by growing ellipses by clouds of points until a
%maximum size is reached, the code will then record this data for each
%point in a solution matrix. It can also establish the maximum radius for
%porosity, and any other parameter which requires a maximum space solution.
%
%Coded by Chris Harbord, RML Durham, last edit: 05/12/2015
%% Some initial parameters
% sigmac=6e9;
% P=normal;
l=size(data1); %extracts the matrix dimensions of the inputted data
%res.area=l(1)*l(2)*dx^2; %outputs the area to the solution array for checking
%disp(['Area of image being analysed: ' num2str(res.area)])
% t=linspace(0,360,5); %set up the vector t to define points around a central location
% points=0;
% r=0; %initialise the ellipse size
%% Solution searching loop
% X=ones(l(1),l(2));
% Y=X;
% yd=linspace(0,l(2)*dx,l(2));
% xd=linspace(0,l(1)*dx,l(1));
%% Dynamic interface calculation %%
for f=1:50;
    sn = (f*5)*1e6;
    ac = sn/sc;
Z1=data1;
data=data1;
void = [];
aspp = [];
%for t=1:20;
%     for j=1:l(1);
%         Z1new(j,:)=zeros(l(2),1);
%         Z2new(j,:)=Z1new(j,:);
%         Z1new(j,1:1000)=Z1(j,end-999:end);
%         Z1new(j,1001:end)=Z1(j,1:end-1000);
%         Z2new(j,1:end-1000)=Z2(j,1001:end);
%         Z2new(j,end-999:end)=Z2(j,1:1000);
%     end
%     Z1=Z1new;
%     Z2=Z2new;
%    data=Z1+Z2;
    ZZ=data;
    topheight=1-(ac); %Get percentage area contact of sample based on normal stress 
    zz=reshape(data,l(1)*l(2),1);%Reshape data to perform CDF
    [FF,XX]=ecdf(zz); %Run Kaplan-Maeir survivor CDF estimate
    n=1; %Initialise counter for asperity extraction 
    for ii=1:numel(FF);
        if FF(ii)>=topheight; %Extract points in top percentage contact 
            cutt(n)=XX(ii); %assign values of height greater than asperity contact height
            n=n+1; %continue counter
        else
        end
    end
    cutt=min(cutt); %Extract minimum of apserity heights in distribution
    data(data>cutt) = 2; %Assign 2's to asperity start
    data(data<=cutt) = 1; %Assign 1's to void start
    data(1)=2; %Assign boundaries to distribution, to eliminate edge effects
    data(end)=2; %Assign boundary
    z=data; %Copy data to new variable to z for processing
%     for i=1:numel(ZZ);
%     if ZZ(i)>cutt;
%         Z1(i)=Z1(i)-((Z1(i)-cutt)/8);
%         Z2(i)=Z2(i)-((Z2(i)-cutt)/8);
%     else
%     end
%     end
    %mapp=xx.*map;
        for j=1:l(2);
            clear map mapx asp
            k=2;
            map=diff(z(:,j));
                if map(1)==-1;
                    mapx=1;
                else
                    mapx=-1;
                end
            asp=1;
            xx=linspace(2,length(map)+1,length(map));
            for i=1:length(xx);
                if map(i)~=0; %identify void/asperity markers
                    asp(k)=xx(i); %assign positions of markers
                    mapx(k)=map(i); %Store asperity/void marker
                    k=k+1;
                else
                end
            end
            if mapx(end)==-1; %Add in boundary based on last marker
            mapx(end+1)=1;
            else
            mapx(end+1)=-1;
            end
            asp(end+1)=xx(end)+1;
            map=diff(asp).*mapx(1:end-1);
            mapv=map;
            mapa=map;
            mapv(mapv>=-1)=NaN;
            mapa(mapa<=1)=NaN;
            mapv(isnan(mapv))=[];
            mapa(isnan(mapa))=[];
            void=horzcat(void, mapv.*dx);
            aspp=horzcat(aspp, mapa.*dx);
%            T(t)=t;
        end
%void=-mean(mean(void));
aspp(aspp<=2e-6)=[]; %Cut results below nyquist
aspp(aspp>=100e-6)=[]; %Cut results above 100 microns

%% Age calculations based on Beeler et al. 2008
% Adding in lines of code to calculate distributions of ages and possible
% temperatures at asperity tips/D0 length calculations.
Tb = 1150; %Melting temperature of mineral
Tf = 150; %Ambient fault temperature
alpha = 7.25e-7; % Thermal diffusivity
pc = 1.97e6; %Thermal capacity
d = [0:1e-6:100e-6]; %Range of sliding velocities
S = 5e9; %Indentation hardness of mineral
V0 = ((pi*alpha)./d).*((pc*(Tb-Tf))/S)^2; %Critical asperity dimensions for melting


%Contact lifetime
V = 1e-6; %Velocity of shearing interface
age = aspp./(2.*V); %State calculation
age0 = mean(age);
Vw = ((pi*alpha)./mean(aspp)).*((pc*(Tb-Tf))/S)^2;
Vweak(f) = Vw;
theta(f) = age0;
normal(f) = sn;
Dc(f) = (mean(aspp)*1e6)/2;
% disp(['D_c prediction= ' num2str(mean(aspp)*1e6) ' microns']);
% disp(['Critical weakening velocity= ' num2str(Vw) 'm/s']); 
% 
% figure(1)
% histogram(aspp.*1e6,[0:1:100],'Normalization','probability');
% xlabel('Asperity length [\mum]');
% ylabel('Count/critical weakening velocity [ms^{-1}]');
% hold on
% plot(d*1e6,V0);
% xlim([0 100])
% 
% figure(2)
% histogram(age,[0:1:100],'Normalization','probability');
% xlabel('\theta [seconds]');
% ylabel('N(\theta)');
% title(['\theta of contacts sliding at ' num2str(V*1e6) '\mum s^{-1}'])
end    
ans.sn=normal;
ans.Vw=Vweak;
ans.Dc=Dc;
ans.theta = theta;

figure(1)
plot(normal./1e6,Dc);
xlabel('Normal stress [MPa]');
ylabel('D_c [\mum]');
 end