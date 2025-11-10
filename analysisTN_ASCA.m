%% Leading example of ANOVA Simultaneous Component Analysis. A
% Tutorial Review. Submitted to Journal of Chemometrics. 2025. 
%
% Please, reference the following paper if using this data:
%
% C. Diaz, C. González-Olmedo, L. Díaz-Beltrán, J. Camacho, P. Mena García, 
% A. Martín-Blázquez, M. Fernández-Navarro, A. L. Ortega-Granados, F. 
% Gálvez-Montosa, J. A. Marchal, et al., “Predicting dynamic response to 
% neoadjuvant chemotherapy in breast cancer: a novel metabolomics 
% approach,” Molecular Oncology, vol. 16, no. 14, pp. 2658–2671, 2022.
%
% Dependencies: 
%
%   - MEDA Toolbox v1.10 at https://github.com/codaslab/MEDA-Toolbox    
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 10/Nov/2025
%
% Copyright (C) 2025  University of Granada, Granada
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Inicialization

clear
close all
clc

data = importdata('TN_data.xlsx');

X = data.data;
obs_l = data.textdata(2:end,1);
class = data.textdata(2:end,2); 
time = X(2:end,1);
tl = {'Basal','Presurgery','Postsurgery'};
timel = tl(time)';
for i=1:length(X(1,:))
    var_l{i} = num2str(X(1,i));
end
X = X(2:end,2:end);
var_l(1) = [];
pac_l = obs_l;
for i=1:length(pac_l)
    pac_l{i} = pac_l{i}(1:5);
end

F = [string(class), string(timel), string(pac_l)];

pvalue = 0.05


%% Check outliers with PCA

[Dst,Qst] = mspcPca(preprocess2D(X,'Preprocess',1),'Preprocessing',2,'PCs',1:2,'ObsClass',class,'ObsLabel',pac_l);
saveas(gcf,'./Figures/pca-mspc')
saveas(gcf,'./Figures/pca-mspc.eps','epsc')


%% Check outliers with ASCA
 
[table, struct] = parglm(preprocess2D(X,'Preprocess',1), F,'Preprocessing',2,'Model', [1 2], 'Fmtc', 0, 'Nested', [1 3], ...
    'Random', [0 0 1], 'Stable', true);
table.Source(2:4)={"Responder","Time","Patient"}

table2latex(table, 'ASCA.tex', '%.2f');

model = asca(struct);

% Anomalies check
q = sum(model.residuals.^2,2);
for f=1:2
    uf = unique(F(:,f));
    g = string;
    for l=1:length(uf)
        g(find(ismember(F(:,f),uf(l)))) = uf(l);
    end
    figure, boxplot(q,g);

    set(gca,'XTickLabel',uf,'XTickLabelRotation',45,'FontSize',18)
    ylabel('Q-statistic of Residuals','FontSize',18)
    saveas(gcf,sprintf('./Figures/Res%d',f))
    saveas(gcf,sprintf('./Figures/Res%d.eps',f),'epsc')
end

f = 3;

uf = unique(F(:,f));
g = zeros(size(q));
for l=1:length(uf)
    g(find(ismember(F(:,f),uf(l)))) = l;
end

[~, ind ] = sort(time);
plotScatter([g(ind),q(ind)],'ObsClass',timel(ind),'EleLabel',pac_l(ind));
set(gca,'XTick',1:length(uf),'XTickLabel',uf,'XTickLabelRotation',45,'FontSize',14)
axis([0 length(uf)+1,0,400])
legend('Basal','Presurgery','Postsurgery','Location','northwest')
ylabel('Q-statistic of Residuals','FontSize',18)
saveas(gcf,sprintf('./Figures/Res%d',f))
saveas(gcf,sprintf('./Figures/Res%d.eps',f),'epsc')

plotScatter([g(ind),q(ind)],'ObsClass',F(ind,1),'EleLabel',pac_l(ind));
set(gca,'XTick',1:length(uf),'XTickLabel',uf,'XTickLabelRotation',45,'FontSize',14)
axis([0 length(uf)+1,0,400])
uf = unique(F(:,1));
legend(uf{1},uf{2},'Location','northwest')
ylabel('Q-statistic of Residuals','FontSize',18)
saveas(gcf,sprintf('./Figures/Res%db',f))
saveas(gcf,sprintf('./Figures/Res%db.eps',f),'epsc')

%% Try with the rank transform
 
[table, struct] = parglm(preprocess2D(rankTransform(X),'Preprocess',1), F,'Preprocessing',2,'Model', [1 2], 'Fmtc', 0, 'Nested', [1 3], ...
    'Random', [0 0 1], 'Stable', true);
table.Source(2:4)={"Responder","Time","Patient"}

table2latex(table, 'ASCA2.tex', '%.2f');

model = asca(struct);

% Anomalies check
q = sum(model.residuals.^2,2);
for f=1:2
    uf = unique(F(:,f));
    g = string;
    for l=1:length(uf)
        g(find(ismember(F(:,f),uf(l)))) = uf(l);
    end
    figure, boxplot(q,g);

    set(gca,'XTickLabel',uf,'XTickLabelRotation',45,'FontSize',18)
    ylabel('Q-statistic of Residuals','FontSize',18)
    saveas(gcf,sprintf('./Figures/Res2%d',f))
    saveas(gcf,sprintf('./Figures/Res2%d.eps',f),'epsc')
end

f = 3;

uf = unique(F(:,f));
g = zeros(size(q));
for l=1:length(uf)
    g(find(ismember(F(:,f),uf(l)))) = l;
end

[~, ind ] = sort(time);
plotScatter([g(ind),q(ind)],'ObsClass',timel(ind),'EleLabel',pac_l(ind));
set(gca,'XTick',1:length(uf),'XTickLabel',uf,'XTickLabelRotation',45,'FontSize',14)
axis([0 length(uf)+1,0,400])
legend('Basal','Presurgery','Postsurgery','Location','northwest')
ylabel('Q-statistic of Residuals','FontSize',18)
saveas(gcf,sprintf('./Figures/Res2%d',f))
saveas(gcf,sprintf('./Figures/Res2%d.eps',f),'epsc')


%% Delete outliers and Check again

indO = find(ismember(F(:,3),'M0357')); % Outlier
X(indO,:) = [];
F(indO,:) = [];
obs_l(indO) = [];
class(indO) = [];
time(indO) = [];
timel(indO) = [];
pac_l(indO) = [];

indO = find(ismember(F(:,3),'M0370')); % Outlier
X(indO,:) = [];
F(indO,:) = [];
obs_l(indO) = [];
class(indO) = [];
time(indO) = [];
timel(indO) = [];
pac_l(indO) = [];

indO = find(ismember(F(:,3),'M0291')); % Outlier
X(indO,:) = [];
F(indO,:) = [];
obs_l(indO) = [];
class(indO) = [];
time(indO) = [];
timel(indO) = [];
pac_l(indO) = [];

[table, struct] = parglm(preprocess2D(X,'Preprocess',1), F,'Preprocessing',2,'Model', [1 2], 'Fmtc', 0, 'Nested', [1 3], ...
    'Random', [0 0 1], 'Stable', true)
table.Source(2:4)={"Responder","Time","Patient"}

table2latex(table, 'ASCA3.tex', '%.2f');

model = asca(struct);

% Anomalies check
q = sum(model.residuals.^2,2);
for f=1:2
    uf = unique(F(:,f));
    g = string;
    for l=1:length(uf)
        g(find(ismember(F(:,f),uf(l)))) = uf(l);
    end
    figure, boxplot(q,g);

    set(gca,'XTickLabel',uf,'XTickLabelRotation',45,'FontSize',18)   
    ylabel('Q-statistic of Residuals','FontSize',18)
    saveas(gcf,sprintf('./Figures/Res3%d',f))
    saveas(gcf,sprintf('./Figures/Res3%d.eps',f),'epsc')
end

f = 3;

uf = unique(F(:,f));
g = zeros(size(q));
for l=1:length(uf)
    g(find(ismember(F(:,f),uf(l)))) = l;
end

[~, ind ] = sort(time);
plotScatter([g(ind),q(ind)],'ObsClass',timel(ind),'EleLabel',pac_l(ind));
set(gca,'XTick',1:length(uf),'XTickLabel',uf,'XTickLabelRotation',45,'FontSize',14)
axis([0 length(uf)+1,0,400])
legend('Basal','Presurgery','Postsurgery','Location','northwest')   
ylabel('Q-statistic of Residuals','FontSize',18)
saveas(gcf,sprintf('./Figures/Res3%d',f))
saveas(gcf,sprintf('./Figures/Res3%d.eps',f),'epsc')


%% Model visualization: only significant terms

% Factor 1

scores(model.factors{1},'ObsLabel',string(pac_l(1:3:end)),'ObsClass',string(class(1:3:end)),'BlurIndex',1);
xlabel('Patients')
saveas(gcf,'./Figures/scoresPCA1')
saveas(gcf,'./Figures/scoresPCA1.eps','epsc')

loadings(model.factors{1},'VarsLabel',var_l,'BlurIndex',1);
xlabel('Metabolomic features')
saveas(gcf,'./Figures/loadingsPCA1')
saveas(gcf,'./Figures/loadingsPCA1.eps','epsc')

% Factor 3
model.factors{3}.lvs=1:2; % Visualize the first 2 PCs
scores(model.factors{3},'ObsLabel',pac_l,'ObsClass',class,'BlurIndex',1);
saveas(gcf,'./Figures/scoresPCA3')
saveas(gcf,'./Figures/scoresPCA3.eps','epsc')

loadings(model.factors{3},'VarsLabel',var_l,'BlurIndex',1);
saveas(gcf,'./Figures/loadingsPCA3')
saveas(gcf,'./Figures/loadingsPCA3.eps','epsc')

% Visualize the MSPC
varPca(model.factors{3}.matrix,'PCs',1:20,'Preprocessing',0); % Either 1 or 4 PCs
mspcPca(model.factors{3}.matrix,'Preprocessing',0, 'PCs',1:4,'PlotCal',false,'ObsTest',model.factors{3}.matrix+model.residuals, 'ObsLabel',pac_l,'ObsClass',class,'PValueD',0.001,'PValueQ',0.000001);
saveas(gcf,'./Figures/residualsPCA3')
saveas(gcf,'./Figures/residualsPCA3.eps','epsc')


% Residuals
modelR = pcaEig(model.residuals,'PCs',1:2); % Visualize the first 2 PCs
scores(modelR,'ObsLabel',pac_l,'ObsClass',class,'BlurIndex',1);
saveas(gcf,'./Figures/scoresRes')
saveas(gcf,'./Figures/scoresRes.eps','epsc')

loadings(modelR,'VarsLabel',var_l,'BlurIndex',1);
saveas(gcf,'./Figures/loadingsRes')
saveas(gcf,'./Figures/loadingsRes.eps','epsc')

% Visualize the MSPC
varPca(model.residuals,'PCs',1:20,'Preprocessing',0);  % 1 PC
mspcPca(model.residuals,'Preprocessing',0, 'PCs',1:1, 'ObsLabel',pac_l,'ObsClass',class);
saveas(gcf,'./Figures/residualsRes')
saveas(gcf,'./Figures/residualsRes.eps','epsc')

%%