clc
clear
close all

%% For demonstration only, the results are not included in our paper
%% NOTE: The following extracted device parameters are different from those used in oAH, oDPI and ONCs
ep_0 = 8.85e-14;    %F.cm^-1
ep_r = 3.15;        %Parylene dielectric constant
d = 200e-7;         %The dielectric thickness in cm
ci = ep_0*ep_r/d;   %Dielectric capacitance
W_n = 1000e-6;      %Channel width
W_p = 1000e-6;      %Channel width
l = 100e-6;         %Channel length
mu0_n = 0.02252;    %Mobility incm2/V.s
mu0_p = 0.3201;     %Mobility incm2/V.s
Vt_n = -10.621;     %Threshold voltage
Vt_p = -4.034;      %Threshold voltage
Vss_n =  3.8383;
Vbb_n = 0.23918;
landa_n = -0.015541;
gama_n = -0.21886;
Vss_p =  0.78061;
Vbb_p = 0.3197;
landa_p = 0.006991;
gama_p = -0.7757;
aaa_n = 5.1444;
bbb_n = 0.16406;
ccc_n = -19.7587;
aaa_p = 0.395;
bbb_p = -0.84;
ccc_p = -23.6763;

C = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],...
    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840],[0 0 0]};
color_list=[204 153 0; 177 119 222; 55 173 107; 241 64 64];
color_list = color_list./255;
Marker_list = ['o','s','^'];
Vw = 15;
Vdd_list = [20,21,22]; Vtau = 15; Csyn = 10; % nF

figure(1);
hs = [];
for i = 1:3
    Vdd = Vdd_list(i); %Input = Input_list{i};
    out = sim('oLDI_exc.slx');
    Isyn = out.Isyn;
    Input = out.Input;
    ax = gca();
    yyaxis left;
    hs(i)=plot(Isyn(:,1),Isyn(:,2),'-','LineWidth',2,'color',C{i});
    hold on; 
    if (i==3)
        ylabel('I_{Syn} (A)','FontSize',18,'FontName','Times New Roman'); ylim([0.8e-10,2.2e-10]);
        yyaxis right;
        hi=plot(Input(:,1),Input(:,2),'-','LineWidth',2,'color','k');
        ylabel('PSP (V)','FontSize',18,'FontName','Times New Roman'); ylim([-23.2,23.2]); xlim([0,50]);
    end
    if (i==3)
        xlabel('Time (s)','FontSize',18,'FontName','Times New Roman');
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        ax.FontSize = 16;
    end
end
legend({'V_{DD} = 20 V','V_{DD} = 21 V','V_{DD} = 22 V','PSP'},'FontSize',14,'FontName','Times New Roman');
% legend boxoff;
% box on;
opt = [];
opt.HoldLines = 1;
a = 4;
b = a/((1+sqrt(5))/2);
opt.BoxDim = [a+b, a]; %[width, height]
opt.FontName = 'Times New Roman';
% opt.LegendLoc = 'NorthWest';
opt.FontSize = 18;
% Save? comment the following line if you do not want to save
opt.FileName = 'figs/oLDI_Vdd.eps'; 
% create the plot
setPlotProp5(opt);
% I0p = -W_p*mu0_p*ci*(Vss_p)^(gama_p+2)/(gama_p+2)/l/(Vbb_p)^gama_p*exp((gama_p+2)/Vss_p*(Vt_p))