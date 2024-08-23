addpath("/Users/haily/Documents/projects/social/startup+code")
addpath("/Users/haily/Documents/projects/social/startup+code/GenLouvain-master")
addpath '/Users/haily/Documents/projects/social/startup+code/fcn'
addpath '/Users/haily/Documents/projects/social'

%% get data
soc = readtable("../mat/completeEnviro_042424.csv");
%soc = soc(:,2:end);
soca = single(table2array(soc(:,2:19)));

%measures = soc.Properties.VariableNames;
%measures = measures(3:end);
%soc.subjectkey = string(erase(soc.subjectkey,"_"));
%soc.eventname = string(soc.eventname);

% save only subs who also have fc data
%load("mat/CompleteFC.mat")

%keepheno = ismember(soc.subjectkey,complete_fc.ids.subjectkey);
%soc_complete = soc(keepheno==1,:);
% save
%writetable(soc_complete,'mat/multilayer_soc_data.csv')

% separate measures into distinct matrices by time point
%soc_complete.eventname = string(soc_complete.eventname);
%soc_base = soc_complete(soc_complete.eventname == 'baseline_year_1_arm_1',:);
%soc2 = soc_complete(soc_complete.eventname == '2_year_follow_up_y_arm_1',:);
%soc3 = soc_complete(soc_complete.eventname == '3_year_follow_up_y_arm_1',:);

% save only people with complete data across time
%complete = intersect(soc_base.subjectkey, intersect(soc2.subjectkey, soc3.subjectkey));
%basekeep = ismember(soc_base.subjectkey,complete);
%sb = soc_base(basekeep==1,:);
%keep2 = ismember(soc2.subjectkey,complete);
%s2 = soc2(keep2==1,:);
%keep3 = ismember(soc3.subjectkey,complete);
%s3 = soc2(keep3==1,:);
%soc_enviro = sb;
%soc_enviro = [soc_enviro; s2];
%soc_enviro = [soc_enviro; s3];
%writetable(soc_enviro,'long_enviro_no_z.csv');
%writematrix(table2array(soc_complete(:,1)),'long_ids.csv');

meas_names = {'fc_1','sd_1','si_1','se_1','pm_1','sr_1','fc_2','sd_2','si_2','se_2','pm_2','sr_2','fc_3','sd_3','si_3','se_3','pm_3','sr_3'};

%% format data for clustering
nmeas = 6;
ntime = 3;
%long = zeros(length(complete),nmeas*ntime);
%nsub = size(long,1);
%long(:,1:nmeas) = table2array(sb(:,3:8));
%long(:,nmeas+1:2*nmeas) = table2array(s2(:,3:8));
%long(:,(2*nmeas)+1:3*nmeas) = table2array(s3(:,3:8));
%writematrix(long,'enviro_longitudinal.csv');
%writematrix(sb(:,1),'ids_longitudinal.csv');
%%% mean subtraction in R, reload
%long = readtable("../mat/long_enviro_noz_minus_mean.csv");
%long = table2array(long);

%% load data that substracts mean from longitudinal scores
%long = readtable('long_enviro_minus_mean.csv');
nsub = height(soc);
L = corr(soca');
L = atanh(L);
L(1:(nsub + 1):end) = 0;

% set params
nreps = 100;   % # louvain repetitions
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

writematrix(Ci,'../mat/ci_long_all.csv');
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

cmapjet = fcn_cmapjet();

% with new order first cluster is 1
end1 = sum(Ci(:,2)==1);

%% make summary figs
%measures = {'Family Environment','School Disengagement','School Involvement','School Environment','Parental Monitoring','Social Resilience'};
%m3 = {'Family Environment','School Disengagement','School Involvement','School Environment','Parental Monitoring','Social Resilience','Family Environment','School Disengagement','School Involvement','School Environment','Parental Monitoring','Social Resilience','Family Environment','School Disengagement','School Involvement','School Environment','Parental Monitoring','Social Resilience'};
m3time = {'Family Environment - Time 1','Family Environment - Time 2','Family Environment - Time 3','School Disengagement - Time 1','School Disengagement - Time 2','School Disengagement - Time 3','School Involvement - Time 1','School Involvement - Time 2','School Involvement - Time 3','School Environment - Time 1','School Environment - Time 2','School Environment - Time 3','Parental Monitoring - Time 1','Parental Monitoring - Time 2','Parental Monitoring - Time 3','Social Resilience - Time 1','Social Resilience - Time 2','Social Resilience - Time 3'}
temporder = [1 1+nmeas 1+(nmeas*2) 2 2+nmeas 2+(nmeas*2) 3 3+nmeas 3+(nmeas*2) 4 4+nmeas 4+(nmeas*2) 5 5+nmeas 5+(nmeas*2) 6 6+nmeas 6+(nmeas*2)];
% subs' scores on measures
%long=table2array(long);
imagesc(soca(order,temporder))
%yline(sum(Ci(:,2)==1),'LineWidth',2)
%yline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2),'LineWidth',2)
%yline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2)+sum(Ci(:,2)==3),'LineWidth',2)
%yline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2)+sum(Ci(:,2)==3)+sum(Ci(:,2)==4),'LineWidth',2)
xticks([1:18])
xticklabels(m3time)
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
%yline(sum(Ci(:,2)==1),':','LineWidth',2)
%yline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2),':','LineWidth',2)
%yline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2)+sum(Ci(:,2)==3),':','LineWidth',2)
%yline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2)+sum(Ci(:,2)==3)+sum(Ci(:,2)==4),':','LineWidth',2)
%xline(sum(Ci(:,2)==1),':','LineWidth',2)
%xline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2),':','LineWidth',2)
%xline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2)+sum(Ci(:,2)==3),':','LineWidth',2)
%xline(sum(Ci(:,2)==1)+sum(Ci(:,2)==2)+sum(Ci(:,2)==3)+sum(Ci(:,2)==4),':','LineWidth',2)
xlabel('subjects'); ylabel('subjects');
colorbar;
colormap(cmapjet)

% assign clusters colors
tgtlvl = 3;
Ci_tgt = Ci(:,tgtlvl);
Ci_tgt(Ci_tgt > 0) = fcn_sort_communities(Ci_tgt(Ci_tgt > 0));
I = fcn_dummyvar(Ci_tgt);
Cent = I'*soca;
cols = distinguishable_colors(size(I,2));
pow = 5;
Ci_col = ones(size(Ci,1),size(Ci,2),3);
for i = 1:size(Ci,2)
    J = fcn_dummyvar(Ci(:,i));
    Centj = J'*soca;
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
save('../mat/ci_hierarch_longi.mat','Ci')