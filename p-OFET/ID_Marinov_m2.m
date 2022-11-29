function z = ID_Marinov_m2(mu0,Vgs_G,Vgs_D,Vt,gama,ci,w,L,...
    landa,Vds_G,Vds_D,Vbb,Vs,Vss,aaa,bbb,ccc)

z1 = zeros(size(Vgs_G,1),size(Vds_G,2));
z2 = zeros(size(Vds_D,1),size(Vgs_D,2));
%IDVG curent calculation
for i = 1:size(Vds_G,2)

    for j = 1:size(Vgs_G,1)
        
        if (Vgs_G(j,i)<=1)
        Vse = Vss*log(1+exp(-(Vgs_G(j,i)-Vt-Vs)/Vss));
        Vde = Vss*log(1+exp(-(Vgs_G(j,i)-Vt-Vds_G(i))/Vss));
        zz = (Vde^(gama+2)-Vse^(gama+2))/(gama+2);
        z1(j,i) = ((w/L)*(mu0/Vbb^gama)*ci*zz*(1+landa*abs(Vds_G(i))));
        else
            z1(j,i) = -exp(aaa*exp(Vgs_G(j,i)*bbb)+ccc);
        end

    end
%     Vse = @(x) (Vss*log(1+exp(-(x-Vt-Vs)/Vss)));
%     Vde = @(x) (Vss*log(1+exp(-(x-Vt-Vds_G(i))/Vss)));
%     f1 = @(x) ((w/L)*(mu0/Vbb^gama)*ci.*(((Vss.*log(1+exp(-(x-Vt-Vds_G(i))/Vss))).^(gama+2)-(Vss.*log(1+exp(-(x-Vt-Vs)/Vss))).^(gama+2))./(gama+2)).*(1+landa*abs(Vds_G(i))));
%     f2 = @(x) -10.^(aaa*x.^bbb+ccc);
%     [dist, x] = meshgrid([0.1],Vgs_G(:,i)');
%     foo = blend(f1,f2,3.5,dist);
%     z1(:,i) = foo(Vgs_G(:,i));
end


%IDVD curent calculation

for i = 1:size(Vgs_D,2)
    
    for j = 1:size(Vds_D,1)
        
        if (Vgs_D(i)<=1)
        Vse = Vss*log(1+exp(-(Vgs_D(i)-Vt-Vs)/Vss));
        Vde = Vss*log(1+exp(-(Vgs_D(i)-Vt-Vds_D(j,i))/Vss));
        zz = (Vde^(gama+2)-Vse^(gama+2))/(gama+2);
        z2(j,i) = ((w/L)*(mu0/Vbb^gama)*ci*zz*(1+landa*abs(Vds_D(j,i))));
        else
            z2(j,i) = -exp(aaa*exp(Vgs_D(i)*bbb)+ccc);
        end

    end
end
z = {z1,z2};
end