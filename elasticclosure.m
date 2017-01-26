    function ans = elasticclosure(data,dx,sn)

%%
%This is a function to calculate the closure of voids in a rough elastic
%surface contacting a flat substrate. It uses the closure stress in
%elliptical cracks to find the minimum stress required to close each point
%through z(x), and uses this to create a new asperity. Stress is estimated
%from elasticity.
% Inputs: 
%       data = surface topography z in m's.
%       dx = sample spacing of data.
%       loops = no of loops that the simulator should perform, please note
%       this shuld be greater than 2 as there is a differencing loop.
%
% Ouputs:
%       avaspp = average of asperity lengths at corresponding stress
%       minaspp = minimum asperity length at corresponding stress
%       maxaspp = largest asperity at corresponding stress
%       stress = stress level corresponding to other values
%       void = maximum void length at corresponding stress level
%       avvoid = average void length, corresponding to lambdac in Scholz 1988
%       lc = estimate critical nucleation length at corresponding normal stress
%       eta = the stability ratio for each corresponding normal stress
%
% CODED BY C.HARBORD, DURHAM RML, LAST EDIT: 24/01/2017
%% Calculation of initial geometry 
tic %start counting time, its a precious commodity
%Sc = 1e9; %Indentation hardness of contacting surface
%Sn = 3e6;
z = data; 
l=size(data); %Extract dimensions of input data
G = 60e9; %Shear modulus
nu = 0.25; %Poisson ratio of elastic medium
k = 4*(1-nu); %Elastic term
gamma = G/2; 
xq = [5e6:5e5:sn]; %Normal stress interpolation array
ans.stressint = xq; %Store normal stress array

for kk = 1:l(2)
clear lc avaspp minaspp maxaspp stress void avvoid nasp
%% Problem initialisation
z = data(:,kk); %Extract 1-d profile in array for elastic contact modelling
[zz,ii] = sort(z); % Sort profile into height to identify highest point in topography
Zasp = zz(end); %Identify elevation of first asperity 
Xasp = ii(end)+1; %Identify position of first asperity 
XaspC = 1; %Assign 1 to provide boundary to problem
x = [1:1:l(1)+2]; %Generate index array for mapping asperity locations
XaspC(2:1+length(Xasp)) = Xasp; %Transfer asperity positions to array with boundary conditions applied
XaspC(2+length(Xasp)) = l(1)+2; %Apply final boundary condition of asperity at end of array
h=waitbar(0,['Iteration= ' num2str(0)]); %Wait bar to visualise progress of array
Zforward = z'; %Store translated topography for use in nect loop
stress=0; %Assign initial stress condition
stressj=0; %Assign initial stress condition indicator
j = 1;
%% Loop to simulate loading of interface
while stressj < sn %Run loop until required level of stress reached
    j = j+1; %Start the loop counter
    waitbar(stressj/sn,h,['Progress= ' num2str((stressj/sn)*100,2) '% of profile ' num2str(kk) ' of' num2str(l(2))]) %Update witbar with latest progress
    clear Xnew Lnew dXasp %Clear variables used to ascertain asperity locations
    dXasp = diff(XaspC); %Difference the grid to obtain the spacing between asperities, L
    dXasp(1) = dXasp(1).*2; %Multiply first value in array by 2 to apply boundary condition
    dXasp(end) = dXasp(end).*2; %Multiply last value in array by 2 to apply boundary condition
    Zasp = data(Xasp); %Extract the elevation of the asperity elevations[not required for solution...
    Xnew=[]; %Initialise Xnew, 
    Lnew=[]; %Initialise Lnew
    for i = 1:length(XaspC)-1;  %Run the asperity positioning loop until all asperities investigated
        clear Lvoid Xval
        %% Asperity positioning
        Xval = XaspC(i+1)-x; %Create translation matrix to assign positions from asperities
        Xval(Xval<1) = []; %Remove results less than 1 as outside domain of void
        Xval(Xval>dXasp(i)-1) = []; %Remove results greater than length of void as outside domain
        Xnew = horzcat(Xnew,(Xval-(dXasp(i)/2)),0); %Concatenate results with new lengths through void domains
        %Void lengths for displacement calculation
        Lvoid = ones(size(Xval)).*(dXasp(i)); %Mark length of void to be solved within
        Lnew = horzcat(Lnew,Lvoid,0); %Concatenate void lengths in sub-domains
    end
    %% Data extraction to remove boundaries
    start=1; %The start of the data range
    finish=l(1); %The end of the data range
    Zmax=max(Zforward); %Find maximum height for establishing closure distances
    DZ=Zmax-Zforward;   %Calculate closure distances 
    %%  Find 
    LnewC = Lnew(start:finish).*dx; %Extract void length information array
    XnewC = Xnew(start:finish).*dx; %Extract distance array
    sigmaX = (DZ.*gamma)./sqrt(LnewC.^2-XnewC.^2); %Calculate stress for closure at each point
    sigmaX([XaspC])=nan; %Replace stress values at asperities with NaN so ignored in stress calcs
    [sigma, sigmaA] = min(sigmaX); %Calculate minimum stress required for closure to form new asperity
    XaspC = horzcat(XaspC,sigmaA); %Include new asperity position in array
    XaspC = sort(XaspC); %Re-sort asperity positions so asperities are in order
    void(j) = max(diff(XaspC)).*pi/2; %Apply geometric statistical correction for 1D data set (email Stefan Nielsen for derivation)
    avvoid(j) = mean(diff(XaspC)).*pi/2; %Apply geometric statistical correction for 1D data set (email Stefan Nielsen for derivation)
    Dz = (sigma./gamma).*sqrt(LnewC.^2-XnewC.^2); %Calculate z displacement update 
    Zforward = Zforward+Dz; %Calculate topography update
    contactstress(j) = (j/l(1)).*10e9; %Calculate contact stress from plastic yield criteria
    stress(j)=stress(j-1)+sigma; %Write stress update, which is cumulative due to incremental stepping
    
%% Asperity calculator from maxellipse6.m
    %This is a loop to calculate the length of sequential asperities using
    %a marker scheme. Lengths are solved by differencing markers.
    clear map mapx asp aspp %Clear variables used for calculations
    x2=ones(size(l(1))); %Assign ones as void markers
    x2([XaspC]) = 2; %Assign asperity marker to new variable x2
    
    
    map=diff(x2); %Difference asperity markers to obtain asperity positions (-1 indicates asp-->void, +1 indicates void-->asperity) 
        if map(1)==-1; %To assign boundary condition of asperity non apserity indicator
            mapx=1; 
        else
            mapx=-1;
        end
    asp=1;
    
    xx=linspace(2,length(map)+1,length(map));
    k=2; %Initialise k
    for i=1:length(xx);
        if map(i)~=0; %identify void/asperity markers
            asp(k)=xx(i); %assign positions of markers
            mapx(k)=map(i); %Store asperity/void marker
            k=k+1; %add 1 to k for next iteration
        else %Else not an asperity so ignore
        end
    end
    if mapx(end)==-1; %Add in boundary based on last marker
        mapx(end+1)=1; 
    else
        mapx(end+1)=-1;
    end
    asp(end+1)=xx(end)+1; %Apply end boundary condition
    map=diff(asp).*mapx(1:end-1); %Multiply differenced asperity map by length vector
    mapv=map; %Store mapped variable to new variable mapv for further calculations
    mapa=map; %Store mapped variable to new variable mapa for further calculations
    mapa(mapa<=1)=NaN; % Replace non asperities with NaN to strip from column vector
    mapa(isnan(mapa))=[]; %Replace nan with 
    aspp=mapa.*dx; %dimensionalise lengths my multiplying by dx
    avaspp(j)=mean(aspp); %Extract average length of asperities
    minaspp(j)=min(aspp); %Extract the minimum dimension of asperity
    maxaspp(j)=max(aspp); %Extract the maximum dimension of asperity
    lc(j)=(mean(aspp)/2*pi()/2*gamma)/(stress(j)*0.65); %Calculate a predicted value of nucleation length given the asperity length and assuming a binary stress distribution
    eta = void./(lc*1e6); %calculate a dimensionless stability ratio eta, based on the maximum void length scaling divided by the critical nucleation length
    stressj = stress(j); %Store the current value of normal stress acting on the surface
    nasp(j) = numel(XaspC); %Calculate the number of asperities on the fault
end
close(h)
ans.avaspp{kk}=avaspp;
ans.minaspp{kk}=minaspp;
ans.maxaspp{kk}=maxaspp;
ans.stress{kk}=stress;
ans.void{kk}=void;
ans.avvoid{kk}=avvoid;
ans.lc{kk}=lc;
ans.eta{kk}=eta;
ans.aspint(kk,:) = interp1(stress(5:end),avaspp(5:end),xq);
ans.voidint(kk,:) = interp1(stress(5:end),void(5:end),xq);
ans.nasp(kk,:) = interp1(stress(5:end),nasp(5:end),xq);

% figure(kk*10+1)
%     plot(stress./1e6,void.*1e-3);
%     xlabel('Normal stress [MPa]');
%     ylabel('Void length [mm]');
%     hold on
%     plot(stress./1e6,avvoid.*1e-3);
%     plot(stress./1e6,lc.*1e3);
%     set(gca,'yscal','log');
%     refline([0 4e2]);
%     legend('Maximum void length','Average void length','Nucleation length','Sample Length');
%     xlim([0 200])
%     ylim([1e-2 1e3])
%     
%     grid on
% figure(kk*10+2)
%     plot(stress./1e6,avaspp.*1e6);
%     hold on
%     plot(stress./1e6,minaspp.*1e6);
%     plot(stress./1e6,maxaspp.*1e6);
%     legend('Average asperity length','Minimum asperity length','Maximum asperity length');
%     xlabel('Normal stress [MPa]');
%     ylabel('Asperity length [\mum]');
%     xlim([0 250])
%     ylim([0 20])
%     grid on
% figure(kk*10+3)
%     plot(stress./1e6,eta);
%     ylabel('Dimensionless stability ratio \eta')
%     xlabel('Normal stress [MPa]')
%     refline([0 1])

end
%% Plotting up the results
figure(1)
plot(xq./1e6,mean(ans.aspint.*1e6))  %Plot averaged asperity lengths as a function of normal stress across all profiles
xlabel('Normal stress [MPa]')
ylabel('Asperity length [\mum]')
xlim([5 sn/1e6])
ylim([0 15])

figure(2)
plot(xq./1e6,mean(ans.voidint.*1e-3)) %Plot average void lengths as a function of normal stress across all profiles
xlabel('Normal stress [MPa]')
ylabel('Void length [mm]')
xlim([5 sn/1e6])
ylim([1e-2 1e3])
set(gca,'yscal','log');
refline([0 4e2]);

figure(3)
plot(xq./1e6,mean(ans.nasp)) %Plot average number of asperities as a function of normal stress across all profiles
xlabel('Normal stress [MPa]')
ylabel('Number of asperities')
xlim([5 sn/1e6])
%This is a joke...
disp(['Quick Donald is coming, best stop doing science and make some alternative facts, by the way this run took ' num2str(toc) ' seconds'])
end
