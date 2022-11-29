function z = CostFunction_general(IDm,IDex_G,IDex_D,gs_G,ds_G,gs_D,ds_D,a,b)

sum1 = 0;
sum2 = 0;

%IDVG calculation
for j = 1:ds_G
    for i = 1:gs_G
        
         sum1 = abs(IDex_G(i,j)-IDm{1,1}(i,j)) + sum1;

    end
%       sum3 = 1/i*(sum1/sum2) + sum3;
%         sum1 = sum1+0*abs(IDm(1,j,1)-IDex_G(1,j))...
%             +1*abs(IDm(end,j,1)-IDex_G(end,j,1))+...
%            +0*abs(IDm(round(gs_G/2),j,1)-IDex_G(round(gs_G/2),j,1));
%         sum1 = sum1+1500*abs(IDm{1,1}(1,j)-IDex_G(1,j))...
%             +100*abs(IDm{1,1}(end,j)-IDex_G(end,j))+...
%            +1500*abs(IDm{1,1}(round(gs_G/2),j)-IDex_G(round(gs_G/2),j)); 
        sum1 = sum1+500*abs(IDm{1,1}(1,j)-IDex_G(1,j))...
            +100*abs(IDm{1,1}(end,j)-IDex_G(end,j))+...
           +1*abs(IDm{1,1}(round(gs_G/2),j)-IDex_G(round(gs_G/2),j));
       sum1 = sum1 + 1000*sum(abs(IDm{1,1}(60:70,j)-IDex_G(60:70,j)));
       sum1 = sum1 + 3000*sum(abs(IDm{1,1}(61,j)-IDex_G(61,j)));
end


for j = 1:gs_D
    for i = 1:ds_D
        
         sum2 = abs(IDex_D(i,j)-IDm{1,2}(i,j)) + sum2;
    end
%     sum2 = sum2+0*abs(IDm(1,j,1)-IDex_D(1,j))...
%             +1*abs(IDm(end,j,1)-IDex_D(end,j,1))...
%             +0*abs(IDm(round(ds_D/2),j,1)-IDex_G(round(ds_D/2),j,1));
% sum2 = sum2+1*abs(IDm{1,2}(1,j)-IDex_D(1,j))...
%             +10*abs(IDm{1,2}(end,j)-IDex_D(end,j))...
%             +10*abs(IDm{1,2}(round(ds_D/2),j)-IDex_D(round(ds_D/2),j));
        sum2 = sum2+1*abs(IDm{1,2}(1,j)-IDex_D(1,j))...
            +1*abs(IDm{1,2}(end,j)-IDex_D(end,j))...
            +10*abs(IDm{1,2}(round(ds_D/2),j)-IDex_D(round(ds_D/2),j));
end

z = a*sum1+b*sum2;
end