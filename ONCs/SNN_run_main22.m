clc
clear
close all

C = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840],[0 0 0]};
%%
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
Vdd_soma = 10; Vdd_syn = 10; %V
Iin1 = 200; Iin2 = 200;% nA
Cmem = 4.7; Cfb=1; Cgs = 1;% nF 
Cload = 0.01;
Vw = 0; Vthr = 5; Csyn = 1; Cgd = 0.8; Vtau = 8; % nF
Vw_exc = 0;
Vw_inh = 0; 
Vthr_inh = 10; %5
Vtau_inh = -5;%-5
Css = 1; %nF
rss = 1; %Ohm

out1 = sim('SNN_circuit22.slx');
Vout = out1.Vout_left;
Vmem = out1.Vmem_left;
Iinput1 = out1.Iinput1;
Iinput2 = out1.Iinput2;
Iinput1(:,2) = Iinput1(:,2)*1e9;
Iinput2(:,2) = Iinput2(:,2)*1e9;

figure(1);
ax = gca();
yyaxis left;
plot(Vout(:,1),Vout(:,2)-0*min(Vout(:,2)),'-','LineWidth',2,'color',C{2});
hold on;
plot(Vmem(:,1),Vmem(:,2)-0*min(Vmem(:,2)),'-','LineWidth',2,'color',C{1});
hold off;
xlabel('Time (s)'); ylabel('V_{out} and V_{mem} (V)');
ylim([-0.1,10]); 
xlim([0,50]);
yyaxis right;
plot(Iinput1(:,1),Iinput1(:,2)-0*min(Iinput1(:,2)),'-','LineWidth',3,'color',C{5});
hold on;
plot(Iinput2(:,1),Iinput2(:,2)-0*min(Iinput2(:,2)),'-','LineWidth',2,'color',C{8});
hold off;
ylabel('I_{in1} and I_{in2} (nA)');
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
box on;
legend({'V_{out}','V_{mem}','I_{in1}','I_{in2}'},'FontSize',14,'FontName','Times New Roman');

opt = [];
opt.HoldLines = 1;
a = 4;
b = a/((1+sqrt(5))/2);
opt.BoxDim = [a+b, a/2]; %[width, height]
opt.FontName = 'Times New Roman';
opt.FontSize = 18;
% Save? comment the following line if you do not want to save
opt.FileName = 'figs/circuit2_1.eps'; 
% create the plot
setPlotProp5(opt);


Iin1 = 300; Iin2 = 200;% nA
Vw=0; Vw_exc=0;
out2 = sim('SNN_circuit22.slx');
Vout = out2.Vout_left;
Vmem = out2.Vmem_left;
Iinput1 = out2.Iinput1;
Iinput2 = out2.Iinput2;
Iinput1(:,2) = Iinput1(:,2)*1e9;
Iinput2(:,2) = Iinput2(:,2)*1e9;
figure(2);
ax = gca();
yyaxis left;
plot(Vout(:,1),Vout(:,2)-0*min(Vout(:,2)),'-','LineWidth',2,'color',C{2});
hold on;
plot(Vmem(:,1),Vmem(:,2)-0*min(Vmem(:,2)),'-','LineWidth',2,'color',C{1});
hold off;
xlabel('Time (s)'); ylabel('V_{out} and V_{mem} (V)');
ylim([-0.1,10]); 
xlim([0,50]);
yyaxis right;
plot(Iinput1(:,1),Iinput1(:,2)-0*min(Iinput1(:,2)),'-','LineWidth',3,'color',C{5});
hold on;
plot(Iinput2(:,1),Iinput2(:,2)-0*min(Iinput2(:,2)),'-','LineWidth',2,'color',C{8});
hold off;
ylabel('I_{in1} and I_{in2} (nA)');
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
box on;
legend({'V_{out}','V_{mem}','I_{in1}','I_{in2}'},'FontSize',14,'FontName','Times New Roman');

opt = [];
opt.HoldLines = 1;
a = 4;
b = a/((1+sqrt(5))/2);
opt.BoxDim = [a+b, a/2]; %[width, height]
opt.FontName = 'Times New Roman';
opt.FontSize = 18;
% Save? comment the following line if you do not want to save
opt.FileName = 'figs/circuit2_2.eps'; 
% create the plot
setPlotProp5(opt);



Iin1 = 200; Iin2 = 200;% nA
Vw=0; Vw_exc=2;
out3 = sim('SNN_circuit22.slx');
Vout = out3.Vout_left;
Vmem = out3.Vmem_left;
Iinput1 = out3.Iinput1;
Iinput2 = out3.Iinput2;
Iinput1(:,2) = Iinput1(:,2)*1e9;
Iinput2(:,2) = Iinput2(:,2)*1e9;
figure(3);
ax = gca();
yyaxis left;
plot(Vout(:,1),Vout(:,2)-0*min(Vout(:,2)),'-','LineWidth',2,'color',C{2});
hold on;
plot(Vmem(:,1),Vmem(:,2)-0*min(Vmem(:,2)),'-','LineWidth',2,'color',C{1});
hold off;
xlabel('Time (s)'); ylabel('V_{out} and V_{mem} (V)');
ylim([-0.1,10]); 
xlim([0,50]);
yyaxis right;
plot(Iinput1(:,1),Iinput1(:,2)-0*min(Iinput1(:,2)),'-','LineWidth',3,'color',C{5});
hold on;
plot(Iinput2(:,1),Iinput2(:,2)-0*min(Iinput2(:,2)),'-','LineWidth',2,'color',C{8});
hold off;
ylabel('I_{in1} and I_{in2} (nA)');
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
box on;
legend({'V_{out}','V_{mem}','I_{in1}','I_{in2}'},'FontSize',14,'FontName','Times New Roman');

opt = [];
opt.HoldLines = 1;
a = 4;
b = a/((1+sqrt(5))/2);
opt.BoxDim = [a+b, a/2]; %[width, height]
opt.FontName = 'Times New Roman';
opt.FontSize = 18;
% Save? comment the following line if you do not want to save
opt.FileName = 'figs/circuit2_3.eps'; 
% create the plot
setPlotProp5(opt);
