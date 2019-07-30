%%%%%%%%%%
%name:          cluster_multi_pairs.m
%description:   clusters regions for multi pairs
%author:        Henriette Miko (henriette.miko@mdc-berlin.de)
%date:          July 18, 2019
%%%%%%%%%%


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
%%%EM%%%%
%now use EM to improve parameters
%keep refining parameters: my and sigma for each node, probabilities of c, 
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

