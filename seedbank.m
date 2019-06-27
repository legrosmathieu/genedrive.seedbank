         function [genfreqDD, genfreqDd, genfreqdd, allelefreqD, allelefreqd, timeto95D] = ...
    seedbank (ngen, mu, c, e, v, m, figson)
%% top-level parameters

% number of generations to run
% ngen=100;

%% create some vectors

% age dependent parameters
muDD=[];
muDd=[];
mudd=[];
bDD=[];
bDd=[];
bdd=[];

% adults
XDD=zeros(ngen,1);
XDd=zeros(ngen,1);
Xdd=zeros(ngen,1);

% seed bank
SDD=zeros(ngen); % time, age
SDd=zeros(ngen);
Sdd=zeros(ngen);

%% parameters

% sex-specific drive conversion bias (0 < e < 1)
% 0 is no drive, 1 is perfect drive
eM=e(1);
eF=e(2);


% genotye-specific fecundity
cDD=c(1); 
cDd=c(2); 
cdd=c(3); 

% genotype-specific age-dependent seed mortality
for age=1:ngen
    muDD(age)=mu(1);
    muDd(age)=mu(2);
    mudd(age)=mu(3);
end

% genotype-specific age-dependent seed germination rate
%
% v parameters = baseline yearly germination rate 
%   (i.e. 1/v is the average time in seed bank ignoring age effects)
%
% m parameters = steepness of decline in germination rate by age

vDD=v(1);
vDd=v(2);
vdd=v(3);
mDD=m(1);
mDd=m(2);
mdd=m(3);
for age=1:ngen
    bDD(age)=vDD/(age.^mDD);
    bDd(age)=vDd/(age.^mDd);
    bdd(age)=vdd/(age.^mdd);
end


%% initialize

timeto95D=-1;

% constitute the initial seed bank 
Sdd(1,1)=1;
for age=2:ngen
    Sdd(1,age)=(1-mudd(age-1))*(1-bdd(age-1))*Sdd(1,age-1);
end

% release drive
XDD(1)=0.05;
Xdd(1)=1-XDD(1);


%% main loop

for t=1:ngen-1
    
    % calculate frequencies of seeds age 0 from adult frequencies
    SDD(t+1,1) = cDD*XDD(t)*(XDD(t)+.5*(1+eM)*XDd(t)) ...
        + .5*(1+eF)*cDd*XDd(t)*(XDD(t)+.5*(1+eM)*XDd(t));
    SDd(t+1,1) = cDD*XDD(t)*(.5*(1-eM)*XDd(t)+Xdd(t)) ...
        + .5*(1+eF)*cDd*XDd(t)*(.5*(1-eM)*XDd(t)+Xdd(t)) ...
        + .5*(1-eF)*cDd*XDd(t)*(.5*(1+eM)*XDd(t)+XDD(t)) ...
        + cdd*Xdd(t)*(XDD(t)+.5*(1+eM)*XDd(t)) ;
    Sdd(t+1,1) = .5*(1-eF)*cDd*XDd(t)*(.5*(1-eM)*XDd(t)+Xdd(t)) ...
        + cdd*Xdd(t)*(.5*(1-eM)*XDd(t)+Xdd(t)) ;
    
    % calculate older seed frequencies 
    for age=2:ngen        
        SDD(t+1,age) = (1-muDD(age-1))*(1-bDD(age-1))*SDD(t,age-1);
        SDd(t+1,age) = (1-muDd(age-1))*(1-bDd(age-1))*SDd(t,age-1);
        Sdd(t+1,age) = (1-mudd(age-1))*(1-bdd(age-1))*Sdd(t,age-1);
    end
    
    % tally new adult frequencies from seed germination
    % replaces previous adult population (annual plants)
    
    tDD=0; 
    tDd=0;
    tdd=0;
 
    for age=1:ngen
        tDD=tDD+(bDD(age)*SDD(t+1,age));
        tDd=tDd+(bDd(age)*SDd(t+1,age));
        tdd=tdd+(bdd(age)*Sdd(t+1,age));
    end
    
    % normalize and update as next generation
    Xsum=tDD+tDd+tdd;
    
    XDD(t+1)=tDD/Xsum;
    XDd(t+1)=tDd/Xsum;
    Xdd(t+1)=tdd/Xsum;
    
    if timeto95D<0 && XDD(t+1)+.5*XDd(t+1)>0.95
        timeto95D=t+1;
    end    
    
end

Stot=SDD+SDd+Sdd;

% output vectors

genfreqDD=XDD;
genfreqDd=XDd;
genfreqdd=Xdd;

allelefreqD=XDD+.5*XDd;
allelefreqd=.5*XDd+Xdd;

meanfitness=XDD*cDD+XDd*cDd+Xdd*cdd;

%% plots

if figson

    % figure 1 - time series of genotype and allele frequencies

    figure1=figure;
    hold on;
    box on;
    x=1:ngen;

    subplot 311;

    plot(x, genfreqDD,'-r', x, genfreqDd, '-g', x, genfreqdd, '-b', 'LineWidth',2);

    set(gca, 'FontSize', 16, 'FontName', 'Calibri');
    xlabel('Generation (year)');
    ylabel('Genotypic frequency');
    legend('DD','Dd','dd','Location', 'East');
    legend boxoff;
   
    subplot 312;

    plot(x, allelefreqD ,'-r', x, allelefreqd, '-b', 'LineWidth',2);

    set(gca, 'FontSize', 16, 'FontName', 'Calibri');
    xlabel('Generation (year)');
    ylabel('Allelic frequency');
    legend('D','d','Location', 'East');
    legend boxoff;
    
    subplot 313;
    
    plot(x, meanfitness, '-k', 'LineWidth', 2);
    xlabel('Generation (year)');
    ylabel('Mean population fitness');
    

    % figure 2 - surface plots of genotype frequencies in seed bank

    figure2=figure;
    box on;
    ymax=ceil(max([2/vDD 2/vDd 2/vdd]));
    y=1:ngen;
    colormap hot;

    tSDD=SDD.';
    tSDd=SDd.';
    tSdd=Sdd.';
    tStot=Stot.';
    zmax=max([max(SDD) max(SDd) max(Sdd)]);

    subplot 311;

    contourf(x,y,tSDD./tStot);
    xlabel('Time (year)');
    ylabel('Age');
    ylim([1 ymax]);
    caxis([0 zmax]);
    title('DD - Genotype frequency');
    set(gca,'Ydir','reverse');

    colorbar('Location','EastOutside');

    hold on;

    subplot 312;

    contourf(x,y,tSDd./tStot);
    xlabel('Time (year)');
    ylabel('Age');
    ylim([1 ymax]);
    caxis([0 zmax]);
    title('Dd - Genotype frequency');
    set(gca,'Ydir','reverse');

    colorbar('Location','EastOutside');

    subplot 313;

    contourf(x,y,tSdd./tStot);
    xlabel('Time (year)');
    ylabel('Age');
    ylim([1 ymax]);
    caxis([0 zmax]);
    title('dd - Genotype frequency');
    set(gca,'Ydir','reverse');

    colorbar('Location','EastOutside');

    % figure 3 - surface plots of allele frequencies in seed bank

    figure3=figure;
    box on;
    colormap hot;

    subplot 211;

    contourf(x,y,(tSDD+.5*tSDd)./tStot);
    xlabel('Time (year)');
    ylabel('Age');
    ylim([1 ymax]);
    caxis([0 zmax]);
    title('D - Allele frequency');
    set(gca,'Ydir','reverse');

    colorbar('Location','EastOutside');

    hold on;

    subplot 212;

    contourf(x,y,(.5*tSDd+tSdd)./tStot);
    xlabel('Time (year)');
    ylabel('Age');
    ylim([1 ymax]);
    caxis([0 zmax]);
    title('d - Allele frequency');
    set(gca,'Ydir','reverse');

    colorbar('Location','EastOutside');


    
end

end


    
    