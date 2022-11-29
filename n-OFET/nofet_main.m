%%%PSO optimization for NOFET parameter extraction based on MARINOV model
%%%By Yi Yang, MJ Mirshojaeian Hosseini and W. Kruger
clc
clear
close all

C = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840],[0 0 0]};
%% Load Experimental Data
IDVG = readmatrix('IDVG.txt');
IDVD = readmatrix('IDVD.txt');

step_G = size(IDVG,1)/sum(IDVG(:) == 0);
step_D = size(IDVD,1)/sum(IDVD(:) == 0);
step = step_G;
paraindex = 2;


for i=1:(size(IDVG,1)/step_G)
    for j = 1:step_G
        data_G1(j,i) = IDVG((paraindex-1)*step_G+j,1);   %Vgs
        result_G1(j,i) = IDVG((paraindex-1)*step_G+j,2); %ID
    end
end
data_G1 = flipud(data_G1);
result_G1 = flipud(result_G1);
clear i j step_G

for i=1:size(IDVD,1)/step_D
    for j = 1:step_D
        data_D1(j,i) = IDVD((i-1)*step_D+j,1);   %Vds
        result_D1(j,i) = IDVD((i-1)*step_D+j,2); %ID
    end
end
clear i j step_D


data_G = data_G1(1:end,:);
result_G = result_G1(1:end,:);
data = flipud(data_G); result = flipud(result_G);

data_D = data_D1(:,:);
result_D = result_D1(:,:);
clear IDVD IDVG

%% OFET parameters
ep_0 = 8.85e-14;    %F.cm^-1
ep_r = 3.15;        %Parylene dielectric constant
d = 200e-7;         %The dielectric thickness in cm
ci = ep_0*ep_r/d;   %Dielectric capacitance
W = 1000e-6;        %Channel width
l = 100e-6;         %Channel length
% gamma and threshold votage are non-intrinsic parameters extracted from
% HVG function (see 10.1109/TED.2009.2033309)
% Threshold & Mobility
F = griddedInterpolant(flipud(data(:,paraindex)),flipud(result(:,paraindex)));
fun = @(t) F(t);
for k = 1:step
    q(k) = integral(fun,data(end,paraindex),data(end-k+1,paraindex));
end
q = q';
HVG = q./flipud(result(:,paraindex));
n = round (size(data,1))/200;
n = round(n);
dt_fit = flipud(data(:,paraindex));
p1 = polyfit(dt_fit(193*n:198*n),HVG(193*n:198*n),1);
Threshold_Marinov = -(p1(2)/p1(1))
% Vt = Threshold_Marinov;
gama = 1/p1(1)-2.87 % For OFET with low Ion/Ioff, we use 2.87 instead of 3
% Note: The following gamma and threshold voltage values (manually tuned) 
% are slightly different from ones given by the HVG function, and also 
% guarantee good fitting performance. Other extracted values are possible.
gama = -0.1937;      %Mobility incm2/V.s
Vt = -0.996;         %Threshold voltage

%Parameters for IDVG
Vds_G = [20];       %Drain voltage for IDVG
Vgs_G = data_G;     %Gate voltage for IDVG
ID_ex_G = result_G;
gs_G = size(Vgs_G,1);
ds_G = size(Vds_G,2);

%Parameters for IDVD
Vgs_D = [6 8 10];   %Gate voltage for IDVD
Vds_D = data_D;     %Drain voltage for IDVD
ID_ex_D = result_D;
gs_D = size(Vgs_D,2);
ds_D = size(Vds_D,1);

Vs = 0;
a = 0.5;  %Weight of IDVG errors in the cost function
b = 0.5;  %Weight of IDVD errors in the cost function
zz = 0;   %BestCost counter
aa = 0;
coe = 1e-3*randi([1 100]);
%% Problem Definition

%CostFunction=@(IDG_m,IDG_ex,gs,ds) err(IDG_m,IDG_ex,gs,ds);   Cost Function
%IDG_m: Vector of ID for model consist of Ids for the range of Vgs and Vds
%IDG_ex: Vector of ID from experimental data consist of Ids for the range of Vgs and Vds
%gs: The size of vector of Vgs
%ds: The size of vector of Vds (the number of data points at specific Vgs)
                                                           

nVar=7;                 % Number of Decision Variables
VarSize=[1 nVar];       % Size of Decision Variables Matrix
VssMin= 3.5005;         % Lower Bound of Vss
VssMax= 3.500501;       % Upper Bound of Vss
VbbMin= 0.45741;        % Lower Bound of Vbb
VbbMax= 0.4574101;      % Upper Bound of Vbb
mu0Min = 0.09278;    % Lower Bound of gamma
mu0Max = 0.0927801;      % Upper Bound of gamma
landaMin = 0.003043;    % Lower Bound of lambda
landaMax = 0.00304301;  % Upper Bound of lambda
aaaMin = 4.322;         % Lower Bound of a in Double-Exponential
aaaMax = 4.32201;       % Upper Bound of a in Double-Exponential
bbbMin = 0.1926;        % Lower Bound of b in Double-Exponential (Combined with thermal voltage)
bbbMax = 0.192601;      % Upper Bound of b in Double-Exponential (Combined with thermal voltage)
cccMin = -20.4301;      % Lower Bound of c in Double-Exponential
cccMax = -20.43;        % Upper Bound of c in Double-Exponential
%% PSO Parameters
MaxIt=500;       % Maximum Number of Iterations
nPop=200;        % Population Size (Swarm Size)
% PSO Parameters
% w=1;           % Inertia Weight
% wdamp=.99;     % Inertia Weight Damping Ratio
% c1=1.5;        % Personal Learning Coefficient
% c2=2;          % Global Learning Coefficient

% % Constriction Coefficients
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;          % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=chi*phi1;    % Personal Learning Coefficient
c2=chi*phi2;    % Global Learning Coefficient
kk = zeros(1,MaxIt);


VelVssMax=(VssMax-VssMin);
VelVssMin=-VelVssMax;
VelVbbMax=(VbbMax-VbbMin);
VelVbbMin=-VelVbbMax;
Velmu0Max=(mu0Max-mu0Min);
Velmu0Min=-Velmu0Max;
VellandaMax=(landaMax-landaMin);
VellandaMin=-VellandaMax;
VelaaaMax=(aaaMax-aaaMin);
VelaaaMin=-VelaaaMax;
VelbbbMax=(bbbMax-bbbMin);
VelbbbMin=-VelbbbMax;
VelcccMax=(cccMax-cccMin);
VelcccMin=-VelcccMax;
%% Initialization
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
particle=repmat(empty_particle,nPop,1);
GlobalBest.Cost=inf;
for i=1:nPop
    
    % Initialize Position
    particle(i).Position(1)=unifrnd(VssMin,VssMax,1);
    particle(i).Position(2)=unifrnd(VbbMin,VbbMax,1);
    particle(i).Position(3)=unifrnd(mu0Min,mu0Max,1);
    particle(i).Position(4)=unifrnd(landaMin,landaMax,1);
    particle(i).Position(5)=unifrnd(aaaMin,aaaMax,1);
    particle(i).Position(6)=unifrnd(bbbMin,bbbMax,1);
    particle(i).Position(7)=unifrnd(cccMin,cccMax,1);


    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    Vss = particle(i).Position(1);
    Vbb = particle(i).Position(2);
    landa = particle(i).Position(4);
    mu0 = particle(i).Position(3);
    aaa = particle(i).Position(5);
    bbb = particle(i).Position(6);
    ccc = particle(i).Position(7);
    ID_m = ID_Marinov_m2(mu0,Vgs_G,Vgs_D,Vt,gama,ci,W,l,...
    landa,Vds_G,Vds_D,Vbb,Vs,Vss,aaa,bbb,ccc);
    particle(i).Cost = CostFunction_general...
        (ID_m,ID_ex_G,ID_ex_D,gs_G,ds_G,gs_D,ds_D,a,b);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end
BestCost=zeros(MaxIt,1);
%% PSO Main Loop
for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position)+zz;
        
        % Apply Velocity Limits
        particle(i).Velocity(1) = max(particle(i).Velocity(1),VelVssMin);
        particle(i).Velocity(1) = min(particle(i).Velocity(1),VelVssMax);
        particle(i).Velocity(2) = max(particle(i).Velocity(2),VelVbbMin);
        particle(i).Velocity(2) = min(particle(i).Velocity(2),VelVbbMax);
        particle(i).Velocity(3) = max(particle(i).Velocity(3),Velmu0Min);
        particle(i).Velocity(3) = min(particle(i).Velocity(3),Velmu0Max);
        particle(i).Velocity(4) = max(particle(i).Velocity(4),VellandaMin);
        particle(i).Velocity(4) = min(particle(i).Velocity(4),VellandaMax);
        particle(i).Velocity(5) = max(particle(i).Velocity(5),VelaaaMin);
        particle(i).Velocity(5) = min(particle(i).Velocity(5),VelaaaMax);
        particle(i).Velocity(6) = max(particle(i).Velocity(6),VelbbbMin);
        particle(i).Velocity(6) = min(particle(i).Velocity(6),VelbbbMax);
        particle(i).Velocity(7) = max(particle(i).Velocity(7),VelcccMin);
        particle(i).Velocity(7) = min(particle(i).Velocity(7),VelcccMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
         
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position(1)<VssMin | particle(i).Position(1)>VssMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        IsOutside=(particle(i).Position(2)<VbbMin | particle(i).Position(2)>VbbMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        IsOutside=(particle(i).Position(3)<mu0Min | particle(i).Position(3)>mu0Max);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        IsOutside=(particle(i).Position(4)<landaMin | particle(i).Position(4)>landaMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        IsOutside=(particle(i).Position(5)<aaaMin | particle(i).Position(5)>aaaMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        IsOutside=(particle(i).Position(6)<bbbMin | particle(i).Position(6)>bbbMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        IsOutside=(particle(i).Position(7)<cccMin | particle(i).Position(7)>cccMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position(1) = max(particle(i).Position(1),VssMin);
        particle(i).Position(1) = min(particle(i).Position(1),VssMax);
        particle(i).Position(2) = max(particle(i).Position(2),VbbMin);
        particle(i).Position(2) = min(particle(i).Position(2),VbbMax);
        particle(i).Position(3) = max(particle(i).Position(3),mu0Min);
        particle(i).Position(3) = min(particle(i).Position(3),mu0Max);
        particle(i).Position(4) = max(particle(i).Position(4),landaMin);
        particle(i).Position(4) = min(particle(i).Position(4),landaMax);
        particle(i).Position(5) = max(particle(i).Position(5),aaaMin);
        particle(i).Position(5) = min(particle(i).Position(5),aaaMax);
        particle(i).Position(6) = max(particle(i).Position(6),bbbMin);
        particle(i).Position(6) = min(particle(i).Position(6),bbbMax);
        particle(i).Position(7) = max(particle(i).Position(7),cccMin);
        particle(i).Position(7) = min(particle(i).Position(7),cccMax);

        % Evaluation
        Vss = particle(i).Position(1);
        Vbb = particle(i).Position(2);
        landa = particle(i).Position(4);
        mu0 = particle(i).Position(3);
        aaa = particle(i).Position(5);
        bbb = particle(i).Position(6);
        ccc = particle(i).Position(7);

        ID_m = ID_Marinov_m2(mu0,Vgs_G,Vgs_D,Vt,gama,ci,W,l,...
    landa,Vds_G,Vds_D,Vbb,Vs,Vss,aaa,bbb,ccc);
    particle(i).Cost = CostFunction_general...
        (ID_m,ID_ex_G,ID_ex_D,gs_G,ds_G,gs_D,ds_D,a,b);
       
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
        
    end
    Cost_Counter(it) = GlobalBest.Cost;
    BestCost(it)=GlobalBest.Cost;

     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ]);

     
    w=w*wdamp;
  
        
    kk = diff(Cost_Counter);                
    if  it > 1
        if kk(it-1) < 1e-2
           
            aa = aa+1;

        end
    end
    
    if it>1      
        if kk(it-1)<1e-2 && aa>20
        coe = 1e-5*randi([1 100]);
        zz(1) = coe*unifrnd(VelVssMin,VelVssMax,1);
        zz(2) = coe*unifrnd(VelVbbMin,VelVbbMax,1);
        zz(3) = coe*unifrnd(Velmu0Min,Velmu0Max,1);
        zz(4) = coe*unifrnd(VellandaMin,VellandaMax,1);
        zz(5) = coe*unifrnd(VelaaaMin,VelaaaMax,1);
        zz(6) = coe*unifrnd(VelbbbMin,VelbbbMax,1);
        zz(7) = coe*unifrnd(VelcccMin,VelcccMax,1);

        elseif kk(it-1)>1e-7 && aa>20
        zz = 0;
        aa = 0;
        end
    end
    
    
    
 end
BestSol = GlobalBest;
%% Results
        Vss =  BestSol.Position(1);
        Vbb = BestSol.Position(2);
        landa = BestSol.Position(4);
        mu0 = BestSol.Position(3);
        aaa = BestSol.Position(5);
        bbb = BestSol.Position(6);
        ccc = BestSol.Position(7);
disp(['Vss = ' num2str(Vss) ', Vbb = ' num2str(Vbb) ', landa = '  num2str(landa) ', mu0 = ' num2str(mu0) ...
    ', aaa = ' num2str(aaa) ', bbb = ' num2str(bbb) ', ccc = ' num2str(ccc)]);
ID_m = ID_Marinov_m2(mu0,Vgs_G,Vgs_D,Vt,gama,ci,W,l,...
    landa,Vds_G,Vds_D,Vbb,Vs,Vss,aaa,bbb,ccc);

figure(1)
for i=1:1:1
   semilogy(data_G(1:end,i),abs(ID_m{1,1}(1:gs_G,i)),'-','color',C{1},'LineWidth',2);
   hold on;
   semilogy(data_G1(1:5:end,i),abs(result_G1(1:5:end,i)),'o','color',C{2},'MarkerSize',10);
end

ylabel('I_{D,n} (A)','FontSize',18,'FontName','Times New Roman');
xlabel('V_{GS} (V)','FontSize',18,'FontName','Times New Roman')
legend('Simulated','Measured','Location','northwest','LineWidth',0.5);


opt = [];
opt.HoldLines = 1;
opt.XLim = [-20, 20]; % [min, max]
a = 4;
b = a/((1+sqrt(5))/2);
opt.BoxDim = [a, a]; %[width, height]
opt.FontName = 'Times New Roman';
opt.FontSize = 18;

% Save? comment the following line if you do not want to save
% supported format: 'eps', 'jpg', 'pdf'
opt.FileName = 'figs/IDVG_NOFET.eps'; 

% create the plot
setPlotProp5(opt);

%% IDVD plot
Markers = ['o','s','^','v','>'];
figure (2)
for i=1:1:gs_D
   plot(data_D(1:end,i),(ID_m{1,2}(1:ds_D,i)),'-','LineWidth',2,'color',C{i});
   hold on;
   plot(data_D(1:5:end,i),ID_ex_D(1:5:end,i),Markers(i),'MarkerSize',10,'color',C{i});
end
ylabel('I_{D,n} (A)','FontSize',18,'FontName','Times New Roman');
xlabel('V_{DS} (V)','FontSize',18,'FontName','Times New Roman')
legend('Simulated','Measured','Location','northwest','LineWidth',0.5);
opt1 = [];
opt1.HoldLines = 1;
opt1.XLim = [-0.1, 20]; % [min, max]
opt1.YLim = [0, 6e-7]; % [min, max]
a = 4;
b = a/((1+sqrt(5))/2);
opt1.BoxDim = [a, a]; %[width, height]
opt1.FontName = 'Times New Roman';
opt1.LegendLoc = 'NorthWest';
opt1.FontSize = 18;

% Save? comment the following line if you do not want to save
opt1.FileName = 'figs/IDVD_NOFET.eps'; 

% create the plot
setPlotProp5(opt1);
%%
figure(3);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
% I0n = W*mu0*ci*(Vss)^(gama+2)/(gama+2)/l/(Vbb)^gama*exp(-(gama+2)/Vss*(Vt))