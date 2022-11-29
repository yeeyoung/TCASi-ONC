clc
clear
close all

%% Load Experimental Data
data = load('Data_oDPI.mat');
data1 = data.EXPVtau8Vw10Vth0Csyn1_;
data2 = data.EXPVtau9Vw10Vth0Csyn1_;
data3 = data.EXPVtau10Vw10Vth0Csyn1_;
Exp1 = {data1,data2,data3};
Input_list = {[data1(:,1), data1(:,3)],[data2(:,1), data2(:,3)],[data3(:,1), data3(:,3)]};

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
C = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840],[0 0 0]};
color_list=[204 153 0; 177 119 222; 55 173 107; 241 64 64];
color_list = color_list./255;
Marker_list = ['o','s','^'];
% oDPI circuit parameters
Vw = 10; Vdd = 10; Vtau_list = [8,9,10]; Vthr = 0; 
Csyn = 1; Cgd_list = [1,0.8,0.8]; Cgs = 2;  % nF

fig=figure(1);
a = 4;
b = a/((1+sqrt(5))/2);
set(gcf,'Position',[300 300 72*(a+b) a/1.3*72]);
hs = [];
hm = [];
hi = [];

time_sft = [6.651+0.358,26.8-5.049+4.258,8.3-5.049+0.558]-0.25+0.1;

ax1 = axes();
 
for i = 1:3
    Cgd = Cgd_list(i);
    Vtau = Vtau_list(i); %Input = Input_list{i};
    out = sim('oDPI_exc.slx');
    Isyn = out.Isyn;
    Input = out.Input;
    hs = plot3(ax1,Isyn(:,1),Vtau*ones(size(Isyn(:,1))),Isyn(:,2),'-','LineWidth',2,'color',C{2});
    view(3);
    hold on; 
    he = plot3(ax1,Exp1{i}(1+time_sft(i)*1000:200:end,1)-time_sft(i),Vtau*ones(size(Exp1{i}(1+time_sft(i)*1000:200:end,1)-time_sft(i))),Exp1{i}(1+time_sft(i)*1000:200:end,2),Marker_list(1),'color',C{1},'MarkerSize',10);
    exper=Exp1{i}(1+time_sft(i)*1000:end,2);
    error = immse(Isyn(1:10000,2),exper(1:10000));
    disp(['error = ',num2str(error)]);
    
end
zlabel(ax1,'I_{Syn} (A)','FontSize',18,'FontName','Times New Roman'); %ylim([0.8e-10,2.2e-10]);
xlabel(ax1,'Time (s)','FontSize',18,'FontName','Times New Roman');
ylabel(ax1,'V_{\tau} (V)','FontSize',18,'FontName','Times New Roman');
xlim(ax1,[0,30]); ylim(ax1,[8,10]); zlim(ax1,[0,4e-7]);
ax2 = axes('Position', ax1.Position);
for i = 1:3
    Vtau = Vtau_list(i);
    hi = plot3(ax2,Input(:,1),Vtau*ones(size(Input(:,1))),Input(:,2),'-','LineWidth',2,'color','k');
    hold on;
end
zlabel(ax2,'Pre-synaptic Signal (V)','FontSize',18,'FontName','Times New Roman'); zlim(ax2,[-20,20]); 
ax1.FontSize = 16;
ax2.FontSize = 16;
z1Lim = zlim(ax1);
z2Lim = zlim(ax2);

linkprop([ax1,ax2], {'XLim','YLim','Position','view'});
% Compute z2-axis ticks to match the position of z1 axis ticks
z1TickNorm = (ax1.ZTick - ax1.ZTick(1) + z1Lim(1)) ./ range(z1Lim);
z2Tick = range(z2Lim)*z1TickNorm + z2Lim(1);
ax2.ZTick = z2Tick; 

ax2.Visible = 'off';
text(ax1,repmat(max(xlim(ax1)),size(ax2.ZTick)), ...
    repmat(min(ylim(ax1)),size(ax2.ZTick)), ...
    ax1.ZTick, strsplit((strtrim(sprintf('%.0f ', ax2.ZTick)))),...
    'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',16);
text(ax1,max(xlim(ax1))*1.17, min(ylim(ax1))*1.00, ...
    mean(z2Lim)+0.4e-7, 'PSP (V)' ,...
    'Rotation', 90,'FontSize',17.6,'FontName','Times New Roman');
% legend({['Sim (V_{w} = ',num2str(Vw_list(1)),' V)'],['Exp (V_{w} = ',num2str(Vw_list(1)),' V)'],...
%     ['Sim (V_{w} = ',num2str(Vw_list(2)),' V)'],['Exp (V_{w} = ',num2str(Vw_list(2)),' V)'],...
%     ['Sim (V_{w} = ',num2str(Vw_list(3)),' V)'],['Exp (V_{w} = ',num2str(Vw_list(3)),' V)'],'Pre-synaptic Signal'},'FontSize',14,'FontName','Times New Roman');
% legend boxoff;
legend(ax2,[hs;he;hi],{'Sim (I_{Syn})','Exp (I_{Syn})','PSP'},'FontSize',14,'FontName','Times New Roman','Location','best','EdgeColor','k','LineWidth',1);
box(ax1,'on');
ax1.LineWidth = 1; ax1.XColor = 'k';ax1.YColor = 'k';ax1.ZColor = 'k';


% I0p = -W_p*mu0_p*ci*(Vss_p)^(gama_p+2)/(gama_p+2)/l/(Vbb_p)^gama_p*exp((gama_p+2)/Vss_p*(Vt_p))