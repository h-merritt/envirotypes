%% check % overlap between baseline and hierarchical clusters

%% load data
% longitudinal clusters, ids, & scores
cil = load("ci_hierarch_longi.mat");
cil=cil.Ci;
pheno_l = readtable("social_6measures_3timepoints_pivot.csv");
% baseline clusters, ids, & scores
load("ci_hierarch.mat");
pheno_b = load('completePheno.mat');
pheno_b = pheno_b.pheno;

%% find subs with data in both
% first longi data includes extra people, so keep only people with baseline
% data, since we know all baseline people have quality imaging data
ids_l = pheno_l(ismember(table2array(pheno_l(:,2)),table2array(pheno_b(:,1))),2);
writetable(ids_l,'ids_long.csv')
ids_b = pheno_b(:,1);
writetable(ids_b,'ids_base.csv')
both = pheno_b(ismember(ids_b,ids_l),:);

%% save community ids of subs in both
ci_b = Ci(ismember(ids_b,ids_l),:);

% compute ari for levels 2 and 3
ari = zeros(2,1);
ari(1,:) = fcn_ari_fast([ci_b(:,2),cil(:,2)]);
ari(2,:) = fcn_ari_fast([ci_b(:,3),cil(:,3)]);
