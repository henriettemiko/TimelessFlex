%%%%%%%%%%
%name:          cluster_pairs_multi_combined.m
%description:   clusters multi regions with combined model for pairs
%author:        Henriette Miko (henriette.miko@mdc-berlin.de)
%date:          July 25, 2019
%%%%%%%%%%


disp(numTimePts);
disp(numTracks)
disp(signalGeneratorDirPromEnh);
disp(numClustersPromEnh),
disp(modelDirPromEnh);
disp(curDir);

%define paths
cd '/fast/AG_Ohler/henriette/PANCREAS_final/bnt-master';
addpath(genpathKPM(pwd));
cd(curDir);


numClusters = numClustersPromEnh + numClustersPromProm + numClustersEnhEnh;
disp(numClusters);

%load models from clustering of multi pairs

PromEnh =  load([modelDirPromEnh, '/', num2str(numClustersPromEnh), '/model-', num2str(numClustersPromEnh), '_afterEM.mat']);
PromProm =  load([modelDirPromProm, '/', num2str(numClustersPromProm), '/model-', num2str(numClustersPromProm), '_afterEM.mat']);
EnhEnh =  load([modelDirEnhEnh, '/', num2str(numClustersEnhEnh), '/model-', num2str(numClustersEnhEnh), '_afterEM.mat']);

disp(PromEnh);
disp(PromProm);
disp(EnhEnh);

%combine 3 models into one combined model


% Prepare DAG structure
% 3 fold changes and 4 histone marks for promoters and for enhancers
% 2*12+1=25 nodes in network
% 25 x 25 matrix (0: no edge, 1: directed edge)

numNodes = (numTimePts*numTracks) + 1;
dag = zeros(numNodes);

% head node has edge to all other nodes
dag(1, 2:numNodes) = 1;

%pancreas trees
% first histone modification tree H3K27ac promoter
dag(2,3) = 1;
dag(3,4) = 1;

% second histone modification tree H3K27me3 promoter
dag(5,6) = 1;
dag(6,7) = 1;

% third histone modification tree H3K4me1 promoter
dag(8,9) = 1;
dag(9,10) = 1;

% 4th histone modification tree H3K4me3 promoter
dag(11,12) = 1;
dag(12,13) = 1;


% first histone modification tree H3K27ac enhancer
dag(14,15) = 1;
dag(15,16) = 1;

% second histone modification tree H3K27me3 enhancer
dag(17,18) = 1;
dag(18,19) = 1;

% third histone modification tree H3K4me1 enhancer
dag(20,21) = 1;
dag(21,22) = 1;

% 4th histone modification tree H3K4me3 enhancer
dag(23,24) = 1;
dag(24,25) = 1;

%%%prepare TAN%%%
disp([datestr(datetime) ': Initializing TAN...']);

discreteNodes = 1; % head node
nodeSizes = ones(1,numNodes);
nodeSizes(1,1) = numClusters;
nodeSizes(1,2:numNodes) = 1;
allnet = mk_bnet(dag, nodeSizes, 'discrete', discreteNodes, 'observed', ... 
[2:numNodes]);


nopre=[ 2 8 14 20 5 11 17 23]

for idx=1:numel(nopre)
   t=nopre(idx)
   
   allmean = [struct(PromEnh.initModel.tan.CPD{t}).mean struct(PromProm.initModel.tan.CPD{t}).mean struct(EnhEnh.initModel.tan.CPD{t}).mean ];
   
allcov=zeros(1,1,numClusters);
allcov(1,1,1:numClustersPromEnh)=struct(PromEnh.initModel.tan.CPD{t}).cov;



allcov(1,1,(numClustersPromEnh+1):(numClustersPromEnh+numClustersPromProm))=struct(PromProm.initModel.tan.CPD{t}).cov;

allcov(1,1,(numClustersPromEnh+numClustersPromProm+1):(numClustersPromEnh+numClustersPromProm+numClustersEnhEnh))=struct(EnhEnh.initModel.tan.CPD{t}).cov;

allweights=double.empty(1,0, numClusters);

allnet.CPD{t} = gaussian_CPD(allnet, t, 'mean', allmean, 'cov', allcov, 'weights', allweights);

end



pre=setdiff(2:numNodes, nopre)

for idx=1:numel(pre)
t=pre(idx)



allmean = [struct(PromEnh.initModel.tan.CPD{t}).mean struct(PromProm.initModel.tan.CPD{t}).mean struct(EnhEnh.initModel.tan.CPD{t}).mean ];
   
allcov=zeros(1,1,numClusters);
allcov=zeros(1,1,numClusters);
   allcov(1,1,1:numClustersPromEnh)=struct(PromEnh.initModel.tan.CPD{t}).cov;
allcov(1,1,(numClustersPromEnh+1):(numClustersPromEnh+numClustersPromProm))=struct(PromProm.initModel.tan.CPD{t}).cov;
allcov(1,1,(numClustersPromEnh+numClustersPromProm+1):(numClustersPromEnh+numClustersPromProm+numClustersEnhEnh))=struct(EnhEnh.initModel.tan.CPD{t}).cov;


allweights=zeros(1,1, numClusters);
   allweights(1,1,1:numClustersPromEnh)=struct(PromEnh.initModel.tan.CPD{t}).weights;
allweights(1,1,(numClustersPromEnh+1):(numClustersPromEnh+numClustersPromProm))=struct(PromProm.initModel.tan.CPD{t}).weights;
allweights(1,1,(numClustersPromEnh+numClustersPromProm+1):(numClustersPromEnh+numClustersPromProm+numClustersEnhEnh))=struct(EnhEnh.initModel.tan.CPD{t}).weights;

allnet.CPD{t} = gaussian_CPD(allnet, t, 'mean', allmean, 'cov', allcov, 'weights', allweights);

end

%init regions
%allnet.CPD{1} = tabular_CPD(allnet, 1, 'CPT', [(3617/9684)*struct(PromEnh.tan.CPD{1}).CPT; (800/9684)*struct(PromProm.tan.CPD{1}).CPT; (5267/9684)*struct(EnhEnh.tan.CPD{1}).CPT]);


%multi regions
allnet.CPD{1} = tabular_CPD(allnet, 1, 'CPT', [(3406/9555)*struct(PromEnh.initModel.tan.CPD{1}).CPT; (687/9555)*struct(PromProm.initModel.tan.CPD{1}).CPT; (5462/9555)*struct(EnhEnh.initModel.tan.CPD{1}).CPT]);


%check if probs sum up to 1
sumprobs=cumsum([(3406/9555)*struct(PromEnh.initModel.tan.CPD{1}).CPT; (687/9555)*struct(PromProm.initModel.tan.CPD{1}).CPT; (5462/9555)*struct(EnhEnh.initModel.tan.CPD{1}).CPT]);
disp(sumprobs);
%%




% head node c is discrete node
% each node is a Gaussian node (except head node), has parameters my and
% sigma
% for each cluster: for each node my and sigma, alpha and beta for each  
% edge (linear regression)



outFile1 = [curDir '/model-' num2str(numClusters) '.mat'];
save(outFile1, 'allnet');



%%%
%now cluster all multi pairs together with combined model
%load multi pairs

promEnh = [signalGeneratorDirPromEnh '/allFold_data.txt'];
disp(promEnh);
promProm = [signalGeneratorDirPromProm '/allFold_data.txt'];
disp(promProm);
enhEnh = [signalGeneratorDirEnhEnh '/allFold_data.txt'];
disp(enhEnh);

pem = importdata(promEnh, '\t');
ppm = importdata(promProm, '\t');
eem = importdata(enhEnh, '\t');
dataOrig = [ pem; ppm; eem ];




%%%Import Data%%%
%dataOrig = importdata(inFile, '\t');
numDataPts = length(dataOrig(:,1));
disp(dataOrig(1,1));
disp(['There are ' num2str(numDataPts) ' data points.']);
disp(dataOrig(1:3,:));
rem = num2cell(dataOrig,1);
rem = cell2num(rem);
data = cell(numDataPts, numNodes);
for (i = 1:numDataPts)
    data(i,2:numNodes) = num2cell(rem(i,:));
end
data = data';
mydata = data;

disp(data(1:25,1:10));
    
%%%%%%%%%%%%%%clustering script


%%%start inference%%%
disp(['Starting Inference...']);
logl=0;
marginal = zeros(numClusters, numDataPts);
for (nd = 1:numDataPts)
    evidence = cell(numNodes, 1);
    evidence(2:numNodes,1) = data(2:numNodes,nd);
    engine = jtree_inf_engine(allnet);
    [engine, ll] = enter_evidence(engine, evidence');
    logl = logl+ll;
    marg = marginal_nodes(engine, 1);
    [none, index] = max(marg.T);
    marginal(:,nd) = marg.T;
    data(1,nd) = num2cell(index); 
end


%%%cluster information%%%
disp(['Writing Results...']);
data = data';
classes = cell2num(data(:,1));
marginal = marginal';
outFile = [curDir '/classes-' num2str(numClusters) '_combinedmodel.txt'];
classe = [classes, marginal];
display(size(classe));
dlmwrite(outFile, classe, '\t');




%%%EM%%%%
%now use EM to improve parameters

% keep refining parameters: my and sigma for each node, probabilities of c, 
%edges are linear regression (alpha, beta)


engine = jtree_inf_engine(allnet);
disp(['Started EM algorithm']);

[allnet, LLtrace, engine] = learn_params_em(engine, mydata, 200, 0.0002);
outFile2 = [curDir '/model-' num2str(numClusters) '_afterEM.mat'];
save(outFile2, 'allnet');

numPar = 0;
for (i = 1:numNodes)
    temp = struct(allnet.CPD{i});
    numPar = numPar + temp.nparams;
end

disp('Model is trained');


%%%%%%%%%%%%%%clustering script

%%%start inference%%%
disp(['Starting Inference...']);
logl=0;
marginal = zeros(numClusters, numDataPts);
for (nd = 1:numDataPts)
    evidence = cell(numNodes, 1);
    evidence(2:numNodes,1) = mydata(2:numNodes,nd);
    engine = jtree_inf_engine(allnet);
    [engine, ll] = enter_evidence(engine, evidence');
    logl = logl+ll;
    marg = marginal_nodes(engine, 1);
    [none, index] = max(marg.T);
    marginal(:,nd) = marg.T;
    mydata(1,nd) = num2cell(index); 
end


%%%cluster information%%%
disp(['Writing Results...']);
mydata = mydata';
classes = cell2num(mydata(:,1));
marginal = marginal';
outFile3 = [curDir '/classes-' num2str(numClusters) '_afterEM.txt'];
classe = [classes, marginal];
display(size(classe));
dlmwrite(outFile3, classe, '\t');


exit

