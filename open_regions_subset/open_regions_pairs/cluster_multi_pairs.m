%%%%%%%%%%
%name:          cluster_multi_pairs.m
%description:   clusters regions for multi pairs
%author:        Henriette Miko (henriette.miko@mdc-berlin.de)
%date:          July 18, 2019
%%%%%%%%%%


disp('here')

disp(numClusters);

disp(numTimePts);
disp(numTracks)
disp(signalGeneratorDir);
disp(numClusters),
disp(modelDir);
disp(curDir);

%define paths
cd '/fast/AG_Ohler/henriette/PANCREAS_final/bnt-master';
addpath(genpathKPM(pwd));
cd(curDir);


%cluster multi pairs with model for init pairs

modelFile=[modelDir, '/2_30/model-', num2str(numClusters), '.mat'];
disp(modelFile);
initModel = load(modelFile);
inFile = [signalGeneratorDir '/allFold_data.txt'];
disp(inFile);

numNodes = (numTimePts*numTracks) + 1;
%%%Import Data%%%
dataOrig = importdata(inFile, '\t');
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
    engine = jtree_inf_engine(initModel.tan);
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
outFile = [curDir '/classes-' num2str(numClusters) '_initmodel.txt'];
classe = [classes, marginal];
display(size(classe));
dlmwrite(outFile, classe, '\t');



%%%%%%%%%%
%now use EM to improve parameters

%% Prepare DAG structure
%% 3 fold changes and 4 histone marks for promoters and for enhancers
%% 2*12+1=25 nodes in network
%% 25 x 25 matrix (0: no edge, 1: directed edge)
%
%numNodes = (numTimePts*numTracks) + 1;
%dag = zeros(numNodes);
%
%% head node has edge to all other nodes
%dag(1, 2:numNodes) = 1;
%
%%pancreas trees
%% first histone modification tree H3K27ac promoter
%dag(2,3) = 1;
%dag(3,4) = 1;
%
%% second histone modification tree H3K27me3 promoter
%dag(5,6) = 1;
%dag(6,7) = 1;
%
%% third histone modification tree H3K4me1 promoter
%dag(8,9) = 1;
%dag(9,10) = 1;
%
%% 4th histone modification tree H3K4me3 promoter
%dag(11,12) = 1;
%dag(12,13) = 1;
%
%
%% first histone modification tree H3K27ac enhancer
%dag(14,15) = 1;
%dag(15,16) = 1;
%
%% second histone modification tree H3K27me3 enhancer
%dag(17,18) = 1;
%dag(18,19) = 1;
%
%% third histone modification tree H3K4me1 enhancer
%dag(20,21) = 1;
%dag(21,22) = 1;
%
%% 4th histone modification tree H3K4me3 enhancer
%dag(23,24) = 1;
%dag(24,25) = 1;
%
%%%%prepare TAN%%%
%disp([datestr(datetime) ': Initializing TAN...']);
%
%discreteNodes = 1; % head node
%nodeSizes = ones(1,numNodes);
%nodeSizes(1,1) = numClusters;
%nodeSizes(1,2:numNodes) = 1;
%tan = mk_bnet(dag, nodeSizes, 'discrete', discreteNodes, 'observed', ... 
%[2:numNodes]);
%tan.CPD{1} = tabular_CPD(tan, 1, 'CPT', 'unif', 'dirichlet_weight', 1, ...
%'dirichlet_type', 'unif');
%for (t = 2:numNodes)
%    tan.CPD{t} = gaussian_CPD(tan, t);
%
%end
%
%%disp(struct(tan.CPD{1}).CPT);
%
%% head node c is discrete node
%% each node is a Gaussian node (except head node), has parameters my and
%% sigma
%% for each cluster: for each node my and sigma, alpha and beta for each  
%% edge (linear regression)
%
%
%%%initialize with cluster assignments and not k-means
%
%date = mydata;
%date(1,:) = num2cell(classes);
%
%%date is fully observed with values for hidden nodes initialized with k-means
%%data is not fully observed, values for hidden nodes are []
%
%tan = learn_params(tan, date);
%%tanFull = tan;
%
%%outFileFull = [curDir '/model-' num2str(numClusters) '_fullbeforeEM.mat'];
%%%save(outFileFull, 'tan');
%disp(['learn_params done']);
%



%%%EM%%%%
% keep refining parameters: my and sigma for each node, probabilities of c, 
%edges are linear regression (alpha, beta)
%data has now cluster assignments
%mydata has no assignments

disp(data(1:25,1:10));
disp(mydata(1:25,1:10));


disp(['Started EM algorithm']);
engine = jtree_inf_engine(initModel.tan);

[initModel.tan, LLtrace, engine] = learn_params_em(engine, mydata, 200, 0.0002);
outFile2 = [curDir '/model-' num2str(numClusters) '_afterEM.mat'];
save(outFile2, 'initModel');

numPar = 0;
for (i = 1:numNodes)
    temp = struct(initModel.tan.CPD{i});
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
    engine = jtree_inf_engine(initModel.tan);
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
outFile = [curDir '/classes-' num2str(numClusters) '_afterEM.txt'];
classe = [classes, marginal];
display(size(classe));
dlmwrite(outFile, classe, '\t');


exit

