% Script for bootstrapping genetic correlations


cd '/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/wtHTT_development/developmental_genes/roi'
%Load files - developmental genes
%first set of genes obtained from the web applications
%dev_genes= { 'TSHZ3','SATB2','TM6SF1','ADRA2A','RORB', 'GPR22','OSTN', 'CYP26A1','LIMCH1','EPHB6','ARHGAP20','MCHR2','KCNQ5','RBM24','FBXW7','SSX2IP','PART1','DOK5','FZD4','SLC26A4'};
%genes below (hashed out) are from the paper
dev_genes = {'CYP26A1','TM6SF1','OSTN','TSHZ3','SSX2IP','FBXW7','SATB2','PART1','MCHR2','RORB'}; 

%Load files - AHBA (structural desikan + cerebellum)
fid2 = 'desikan_fsl_cerebellum_exp.csv';
importdata(fid2);
AHBA_ind=ans.data(:,1);
%AHBA_data=ans.data(1:110,2:end);
    %trying only with cortical ROI
AHBA_data=ans.data(1:110,2:end);
AHBA_genes=ans.textdata';
AHBA_genes(1,:) = [];
clear ans

% AHBA  - removes NaNs and centre with z-score
m.nan_rows = find(all(isnan(AHBA_data),2)); 
AHBA_data_nonan_zscore = zscore(AHBA_data(all(~isnan(AHBA_data),2),:));

% Match data 
match_idx = find(ismember(dev_genes,AHBA_genes)); % index refers to Allen data index
match_data = AHBA_data_nonan_zscore(:,match_idx);
HTT_match_idx = find(strcmp(AHBA_genes, 'HTT'));
HTT_match_data = AHBA_data_nonan_zscore(:,HTT_match_idx);

% % PCA on  dev gene lists
[~,PC1] = pca(match_data);

[coeff,score,latent,tsquared,explained,mu] = pca(match_data);
Xcentered = score*coeff'
biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'CYP26A1','TM6SF1','OSTN','TSHZ3','SSX2IP','FBXW7','SATB2','PART1','MCHR2','RORB'});