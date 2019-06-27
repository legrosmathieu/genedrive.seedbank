ngen=200;
mu=[0.1 0.1 0.1];
c=[1.0 1.0 1.0];
%e=[0.9 0.9];
%v=[0.2 0.2 0.2];
m=[0.0 0.0 0.0];
figson=false;

vv=0.1:0.02:0.9;
ee=0.02:0.02:1;

tt=zeros(length(vv),length(ee));
afD=zeros(length(vv),length(ee));

for v=1:length(vv)
    for e=1:length(ee)
        
        [~, ~, ~, allelefreqD, ~, timeto95D] = ... 
            seedbank(ngen, mu, c, [ee(e) ee(e) ee(e)], [vv(v) vv(v) vv(v)], m, figson);
        
        if timeto95D<0
            tt(v,e)=ngen;
        else
            tt(v,e)=timeto95D;
        end
        afD(v,e)=allelefreqD(ngen);

    end
end


%% plot

figure1=figure;
set(figure1, 'Position', [0 0 800 800])

box on;
colormap hot;

contourf(ee,1./vv,tt);
ylabel('1/v - Expected time spent in seed bank (years)');
xlabel('e - Drive conversion rate');

cb=colorbar;

set(gca, 'FontSize', 16, 'FontName', 'Calibri');
set(get(cb, 'Label'), 'String', 'Time to 95% drive allele');
