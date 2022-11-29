clc
clear
close all

C = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840],[0 0 0]};
%% Load Experimental Data
data = load('Data_oAH.mat');
Exp1 = data.ExpVdd10Iin200Cmem47Cfb100_;

ep_0 = 8.85e-14;     %F.cm^-1
ep_r = 3.15;         %Parylene dielectric constant
d = 200e-7;          %The dielectric thickness in cm
ci = ep_0*ep_r/d;    %Dielectric capacitance
W_n = 1000e-6;       %Channel width of n-OFET
W_p = 1000e-6;       %Channel width of p-OFET
l = 100e-6;          %Channel length of p and n-OFET
mu0_n = 0.09278;     %Mobility incm2/V.s of n-OFET
mu0_p = 0.6133;      %Mobility incm2/V.s of p-OFET
Vt_n = -0.996;       %Threshold voltage of n-OFET
Vt_p = -3.5411;      %Threshold voltage of p-OFET
Vss_n =  3.5005;     %Vss of n-OFET
Vbb_n = 0.45741;     %Vbb of n-OFET
landa_n = 0.003043;  %lambda of n-OFET
gama_n = -0.1937;    %gamma of n-OFET
Vss_p =  1.203;      %Vss of p-OFET
Vbb_p = 0.5801;      %Vbb of p-OFET
landa_p = 0.002387;  %lambda of p-OFET
gama_p = -0.7553;    %gamma of p-OFET
aaa_n = 4.32201;     %a_n for double-exponential
bbb_n = 0.1926;      %b_n for double-exponential
ccc_n = -20.4301;    %c_n for double-exponential
aaa_p = 3.1701;      %a_p for double-exponential
bbb_p = -0.5746;     %b_p for double-exponential
ccc_p = -23.0701;    %c_p for double-exponential

color_list=[204 153 0; 177 119 222; 55 173 107; 241 64 64];
color_list = color_list./255;
marker_list = ['-','-','-','-'];
Marker_list = ['o','s','^'];
Vdd = 10; %V
Iin = 200;% nA
Cmem = 47; Cfb=100; Cgs = 0.01;% nF 
Cload = 1e-2;
Css = 1; %nF
rss = 1; %Ohm

out = sim('Organic_AH.slx');
Vout = out.Vout;
Vmem = out.Vmem;

ax = gca();
plot(Vout(:,1),Vout(:,2)-0*min(Vout(:,2)),'-','LineWidth',2,'color',C{2});
hold on;
plot(Exp1(1:100:end,1),Exp1(1:100:end,3),Marker_list(1),'color',C{2},'MarkerSize',10);
exper=Exp1(1:end,3);
error = immse(Vout(12500:25000,2),exper(12500:25000));
disp(['error = ',num2str(error)]);
error1 = error/immse(zeros(12501,1),exper(12500:25000));
disp(['error1 = ',num2str(error1)]);
plot(Vmem(:,1),Vmem(:,2)-0*min(Vmem(:,2)),'-','LineWidth',2,'color',C{1});
plot(Exp1(1:100:end,1),Exp1(1:100:end,2),Marker_list(2),'color',C{1},'MarkerSize',10);
exper=Exp1(1:end,2);
error = immse(Vmem(12500:25000,2),exper(12500:25000));
disp(['error = ',num2str(error)]);
error1 = error/immse(zeros(12501,1),exper(12500:25000));
disp(['error1 = ',num2str(error1)]);
hold off;
xlabel('Time (s)'); ylabel('V_{out} and V_{mem} (V)');
xlim([0,50]);
box on;
legend({'Sim (V_{out})','Exp (V_{out})','Sim (V_{mem})','Exp (V_{mem})'},'FontSize',14,'FontName','Times New Roman');
opt = [];
opt.HoldLines = 1;
a = 4;
b = a/((1+sqrt(5))/2);
opt.BoxDim = [a+b, a/2]; %[width, height]
opt.FontName = 'Times New Roman';
opt.FontSize = 18;
% Save? comment the following line if you do not want to save
opt.FileName = 'figs/oAH_exp4.eps'; 
% create the plot
setPlotProp5(opt);