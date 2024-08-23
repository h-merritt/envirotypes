clear all
close all
clc

%% get hierarchical clusters of social environment data
% rick + haily: spring 2023 - summer 2024
% this script:
% 1. take in pheno data, compute concordance matrix
% 2. perform hierarchical clustering
% 3. create figs of clusters

%% load data, set vars
pheno = readtable('completePheno_9measures.csv'); 
% pheno has baseline scores for all subs with complete data on measures of
% interest
vars = pheno.Properties.VariableNames; % colnames
meas = table2array(pheno(:,[5,7:14])); 
% extract only social environment measures for clustering, save in meas
id_age_sex = pheno(:,[2:4]); % extract id, age, sex
[nsub,nfeat] = size(meas); % get dims for later

% calculate concordance matrix
L = corr(meas');
L = atanh(L);
L(1:(nsub + 1):end) = 0;

%% parameters for clustering
nreps = 1000;   % # louvain repetitions
nrand = 10000;  % # of null models
pcrit = 0.05;   % # statistical criterion
minsize = 1;    % keep clusters that are at least this big
minscan = 1;    % keep clusters that al
G{1} = L;
OrigIndx{1} = (1:length(L))';
goflag = true;

%% loop over hierarchical levels
lvl = 1;
fprintf('level %i, %i community/communities\n',lvl,length(G));
Ci(:,lvl) = ones(length(L),1);

% until termination criterion is reached, keep clustering
while goflag
    
    % update level index
    lvl = lvl + 1;
    
    % set goflag to false - if any sub-communities are detected we'll set
    % it to true later on
    goflag = false;
    
    % for storing labels
    Gnew = [];
    OrigIndxNew = [];
    Samples = [];
    
    % reset count
    count = 0;
    c = zeros(length(L),1);
    fprintf('level %i, %i community/communities\n',lvl,length(G));
    tic;
    
    % loop over all current sub-communities
    for i = 1:length(G)
        %% set up modularity matrix
        C = G{i};
        N = length(C);
        P = nanmean(C(triu(ones(N),1) > 0));
        B = (C - P).*~eye(N);
        %% sample communities
        ci = zeros(N,nreps);
        for irep = 1:nreps
            ci(:,irep) = fcn_community_unif(B);
        end
        Samples{i} = ci; % we don't end up doing anything with this
        %% calculate similarity
        cicon = fcn_consensus_communities_hier(ci,nreps,false,'bct');
        %% calculate modularity contributions
        I = dummyvar(cicon);
        qc = diag(I'*(B*I));
        %% run null model
        qcr = zeros(length(qc),nrand);
        for j = 1:nrand
            r = randperm(N);
            Ir = I(r,:);
            qcr(:,j) = diag(Ir'*(B*Ir));
        end
        %% calculate community size and number of scans contained
        numscans = ones(size(I,2),1);
        commsize = numscans;
        for j = 1:size(I,2)
%             vals = id(I(:,j) > 0);
%             numscans(j) = length(unique(vals));
            commsize(j) = sum(I(:,j));
        end
        %% calculate p-values
        p = mean(bsxfun(@ge,qcr,qc),2);
        idx = find(p < pcrit & numscans >= minscan & commsize >= minsize);
        for j = 1:length(idx)
            count = count + 1;
            jdx = cicon == idx(j);
            Gnew{count} = C(jdx,jdx);
            OrigIndxNew{count} = OrigIndx{i}(jdx);
            c(OrigIndxNew{count}) = count;
        end
        fprintf('%i->%i (total time = %.2f s)\n',i,length(idx),toc);
    end
    if count > 0
        goflag = true;
        G = Gnew;
        OrigIndx = OrigIndxNew;
    end
    Ci(:,lvl) = c;
end
%% organize outputs

I = [];
for j = 1:lvl
    I = [I,fcn_dummyvar(Ci(:,j))];
end
d = I*I';
e = 1 - (d/lvl);
pd = e(tril(ones(length(e)),-1) > 0)';
ll = linkage(pd,'average');
order = optimalleaforder(ll,pd);

save('baseline_ci_all.mat','Ci');


cmapjet = fcn_cmapjet();

%% make a few summary figures
measures = {'Community Risk', 'Neighborhood Safety', 'Parental Monitoring', 
    'Family Conflict', 'Caregiver Acceptance', 'School Engagement', 
    'School Involvement', 'School Disengagement', 'Social Resilience'};
% subs' scores on measures
imagesc(pheno(order,:))
xticks([1:14])
xticklabels(measures)
xlabel('measures')
xtickangle(45)
ylabel('subjects')
colorbar
colormap(cmapjet)
caxis([-3 3])


% coassignment matrix
figure('units','inches','position',[2,2,5,5]);
imagesc(d(order,order))
xlabel('subjects'); ylabel('subjects');
colorbar;

% similarity by clusters
figure('units','inches','position',[2,2,5,5]);
imagesc(L(order,order),[-1,1])
xlabel('subjects'); ylabel('subjects');
colorbar;
colormap(cmapjet)

% assign clusters colors
tgtlvl = 3;
Ci_tgt = Ci(:,tgtlvl);
Ci_tgt(Ci_tgt > 0) = fcn_sort_communities(Ci_tgt(Ci_tgt > 0));
I = fcn_dummyvar(Ci_tgt);
Cent = I'*phenot;
cols = distinguishable_colors(size(I,2));
pow = 5;
Ci_col = ones(size(Ci,1),size(Ci,2),3);
for i = 1:size(Ci,2)
    J = fcn_dummyvar(Ci(:,i));
    Centj = J'*phenot;
    dd = 1 - pdist2(Cent,Centj,'correlation');
    for j = 1:size(J,2)
        w = dd(:,j).^pow;
        w = w.*(w > 0);
        w = w/sum(w);
        e = sum(bsxfun(@times,cols,w));
        Ci_col(J(:,j) > 0,i,:) = repmat(e,sum(J(:,j)),1);
    end
end

% visualize communities hierarchically
figure('units','inches','position',[2,2,5,5]);
imagesc(Ci_col(order,:,:));
xlabel('hierarchical levels');
ylabel('subjects');

%% save
save('../mat/ci_hierarch.mat','Ci')