clearvars;
clc;
close all;

%%DE: Dispersion Estimates; 
%%FE: Focalization Estimates;
%%ASE: Articulatory Space Estimates

%% ===================================================================== %%
load('Tab_Language532Data.mat');
TableData = Tab_Language532Data; 
clear Tab_Language532Data
Region = table2cell(TableData(:,13));
Energy = table2array(TableData(:,10:12)); 
DE = Energy(:, 1); FE = Energy(:, 2);  ASE = Energy(:, 3);
EffectiveDE = DE./ASE;
FE = FE;
logEffectiveDE = log10(DE./ASE); %% -1*log10(DE./ASE); DE./ASE;
logFE = log10(Energy(:, 2));  %% -1*log10(Energy(:, 2)); Energy(:, 2);

%% ===================================================================== %%
%% Fit Power Law model for original valuesof EffectiveDE andn FE
% figure;
% [graphType, slope, intercept, MSE, R2, S] = logfit(EffectiveDE, FE, 'loglog');


figure;
subplot(2,1,1)
gscatter(EffectiveDE, FE, Region);
hold on;
createFit(EffectiveDE, FE)
xlabel('Effective DE'); ylabel('FE');
title('Effective DE VS. FE')
hold off;


subplot(2,1,2);
h = gscatter(logEffectiveDE, logFE, Region);
legend(h, 'Location', 'NorthEast' )
hold on;
mdl = fitlm(logEffectiveDE, logFE);
plot(mdl);
title('log10(Effective DE) VS. log10(FE)')
hold off
xlabel('log10(Effective DE)');
ylabel('log10(FE)');
grid on;

aa=1;  %%For debug


%% ===================================================================== %%
%% Geographic correaltion among regions
tbl = TableData;
Region = table2cell(tbl(:,13));
[Region_Name,ia,ic] = unique(Region);
Region_EffectiveDE = zeros(length(ia),1); 
Region_FE = zeros(length(ia),1);
figure;
for i = 1:length(Region_Name)
    clear tempEnergy tempOptimization temptbl
    idxtemp = (ic==i);
    temptbl = tbl(idxtemp, :);
    Region_Name(i)
    N(i) = sum(idxtemp);
    tempEffectiveDE = EffectiveDE(idxtemp);
    tempFE = FE(idxtemp);
    [rho_DE_FE(i), pval_DE_FE(i)] = corr(tempEffectiveDE, tempFE, 'type', 'Spearman'); %%

    subplot(3,3,i);
    mdl = fitlm(log10(tempEffectiveDE), log10(tempFE));
    h2 = plot(mdl); 
    slope = table2array(mdl.Coefficients(2,1));
    p_value = table2array(mdl.Coefficients(2,4));
    str = strcat(Region_Name(i,1), ', Slope K = ', num2str(slope), ', p-value = ', num2str(p_value));
    title(str);    

    
    Region_EffectiveDE(i) = median(logEffectiveDE(idxtemp));
    Region_EffectiveDE_std(i) = std(logEffectiveDE(idxtemp));
    Region_FE(i) = median(logFE(idxtemp));
    Region_FE_std(i) = std(logFE(idxtemp));
end

aa=1;  %%For debug


%% ===================================================================== %%
%% Correlation between EffectiveDE and FE according to the stratification of language families
tbl = TableData;
clear tempVowelNum tempEffectiveDE tempFE
class = tabulate(table2cell(tbl(:,5)));
[B,I] = sort(cell2mat(class(:,2)),'descend');
class = class(I,:);

num = 12; %% The first number of language families
topNum = table;
rr = zeros(num,1);
pval = zeros(num,1);
figure;
k = 1;
for i=1: num
    class(i,1)
    idxtemp = strcmp(table2cell(tbl(:,5)),class(i,1));
    tempEffectiveDE = logEffectiveDE(idxtemp);
    tempFE = logFE(idxtemp);
    [rr(i),pval(i)] = corr(tempEffectiveDE,...
        tempFE,'type','Spearman');
    ttt(i) = size(tempFE,1);
    combtemp = table(tempEffectiveDE, tempFE);
    topNum = [topNum; [combtemp, tbl(idxtemp, 5)]];
    clear combtemp

    if i~=6  %the sixth rank is the list of isolation languages.
        subplot(4,3,k);
        mdl = fitlm(tempEffectiveDE, tempFE);
        plot(mdl); 
        slope = table2array(mdl.Coefficients(2,1));
        p_value = table2array(mdl.Coefficients(2,4));
        str = strcat(class(i,1), ', Slope K = ', num2str(slope), ', p-value = ', num2str(p_value));
        title(str);

        Family_EffectiveDE(k) = median(logEffectiveDE(idxtemp));
        Family_EffectiveDE_std(k) = std(logEffectiveDE(idxtemp));
        Family_FE(k) = median(logFE(idxtemp));
        Family_FE_std(k) = std(logFE(idxtemp));
        Family_Name(k) = class(i,1);

        k = k+1;
    end
end

aa=1;  %%For debug


%%=======Median values in each groups and plot the results===============%%
figure;
subplot(1,2,1);
[coeff_Region, pval_Region] = corr(Region_EffectiveDE, Region_FE,'type','Spearman');
mdl = fitlm(Region_EffectiveDE, Region_FE);
plot(mdl); 
hold on;
gscatter(Region_EffectiveDE, Region_FE, Region_Name);
hold off;
title('Cross Geographic Regions');
xlabel('log10(Effective DE)');
ylabel('log10(FE)');

subplot(1,2,2);
[coeff_Family, pval_Family] = corr(Family_EffectiveDE, Family_FE,'type','Spearman');
mdl = fitlm(Family_EffectiveDE, Family_FE);
plot(mdl); 
hold on;
gscatter(Family_EffectiveDE', Family_FE', Family_Name');
hold off;
title('Cross Language Families');
xlabel('log10(Effective DE)');
ylabel('log10(FE)');




