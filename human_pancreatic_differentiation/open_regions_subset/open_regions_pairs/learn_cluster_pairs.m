%%%%%%%%%%
%name:          learn_cluster_pairs.m
%description:   learns model and clusters regions for pairs
%author:        Henriette Miko (henriette.miko@mdc-berlin.de)
%date:          July 17, 2019
%%%%%%%%%%


disp(startNum);
disp(endNum);
disp(numTimePts);
disp(numTracks)
disp(signalGeneratorDir);
disp(curDir);


%define paths
cd '/fast/AG_Ohler/henriette/PANCREAS_final/bnt-master';
addpath(genpathKPM(pwd));
cd(curDir);


for (numClusters = startNum:endNum)
    disp(numClusters);

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


    inFileFiltered = [signalGeneratorDir '/allFold_filtered_data.txt'];
    inFile = [signalGeneratorDir '/allFold_data.txt'];
    disp(inFileFiltered);
    disp(inFile);

    disp(numNodes);

    %%%prepare TAN%%%
    disp([datestr(datetime) ': Initializing TAN...']);

    discreteNodes = 1; % head node
    nodeSizes = ones(1,numNodes);
    nodeSizes(1,1) = numClusters;
    nodeSizes(1,2:numNodes) = 1;
    tan = mk_bnet(dag, nodeSizes, 'discrete', discreteNodes, 'observed', ... 
    [2:numNodes]);
    tan.CPD{1} = tabular_CPD(tan, 1, 'CPT', 'unif', 'dirichlet_weight', 1, ...
    'dirichlet_type', 'unif');
    for (t = 2:numNodes)
        tan.CPD{t} = gaussian_CPD(tan, t);

    end

    %disp(struct(tan.CPD{1}).CPT);

    % head node c is discrete node
    % each node is a Gaussian node (except head node), has parameters my and
    % sigma
    % for each cluster: for each node my and sigma, alpha and beta for each  
    % edge (linear regression)

    %%%Import Data%%%
    dataOrigFiltered = importdata(inFileFiltered, '\t');
    numDataPtsFiltered = length(dataOrigFiltered(:,1));
    disp(dataOrigFiltered(1,1));
    disp(['There are ' num2str(numDataPtsFiltered) ' data points.']);
    disp(dataOrigFiltered(1:3,:));
    remFiltered = num2cell(dataOrigFiltered,1);
    remFiltered = cell2num(remFiltered);
    dataFiltered = cell(numDataPtsFiltered, numNodes);
    for (i = 1:numDataPtsFiltered)
        data(i,2:numNodes) = num2cell(remFiltered(i,:));
    end
    dataFiltered = dataFiltered';


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

    disp(data(1:25,1:10));

    %%%Initialize with k-means%%%
    % initial clustering: assign promoter regions to clusters 1,..., 11
    % get initial parameters for EM algorithm

    %(k-1)-means for filtered data

    disp(['Initializing with K-means']);
    k = kmeans(dataOrigFiltered, (numClusters-1), 'Distance', 'cityblock', ...
    'Replicates', 15, 'MaxIter', 300, 'display', 'final');

    koutFile = [curDir '/k-means_initalization_' num2str(numClusters) '.txt'];
    dlmwrite(koutFile, k, '\t');


    numDataPtsNoise = numDataPts - numDataPtsFiltered;
    disp(numDataPtsNoise);

    %noise regions get assigned to cluster k
    NoiseAssign = numClusters+zeros(numDataPtsNoise,1);
    kFull = [k; NoiseAssign];

    disp(length(kFull));


    date = data;
    date(1,:) = num2cell(kFull);

    %date is fully observed with values for hidden nodes initialized with 
    %k-means
    %data is not fully observed, values for hidden nodes are []

    tan = learn_params(tan, date);
    tanFull = tan;

    outFileFull = [curDir '/model-' num2str(numClusters) '_fullbeforeEM.mat'];
    %save(outFileFull, 'tan');
    disp(['learn_params done']);


    %%%EM%%%%

    % keep refining parameters: my and sigma for each node, probabilities of c, 
    %edges are linear regression (alpha, beta)

    engine = jtree_inf_engine(tan);
    disp(['Started EM algorithm']);

    [tan, LLtrace, engine] = learn_params_em(engine, data, 200, 0.0002);
    outFile = [curDir '/model-' num2str(numClusters) '.mat'];
    save(outFile, 'tan');

    numPar = 0;
    for (i = 1:numNodes)
        temp = struct(tan.CPD{i});
        numPar = numPar + temp.nparams;
    end

    disp('Model is trained');


    %store likelihood and number of parameters for likelihood ratio test later
    likelihood = LLtrace(length(LLtrace));
    fname3 = [curDir '/likelihood_' num2str(numClusters) '.mat'];
    f3 = fopen(fname3, 'at');
    cold3 = [num2str(numClusters) '\t' num2str(likelihood, '%.6f')];
    fprintf(f3, cold3);
    fprintf(f3, '\n');
    fclose(f3)

    fname4 = [curDir '/numPar_' num2str(numClusters) '.mat'];
    f4 = fopen(fname4, 'at');
    cold4 = [num2str(numClusters) '\t' num2str(numPar, '%.6f')];
    fprintf(f4, cold4);
    fprintf(f4, '\n');
    fclose(f4)

    AIC = (-2 * LLtrace(length(LLtrace))) + (numPar * 2);
    fname = [curDir '/AIC_' num2str(numClusters) '.mat'];
    f = fopen(fname, 'at');
    cold = [num2str(numClusters) '\t' num2str(AIC, '%.6f')];
    fprintf(f, cold);
    fprintf(f, '\n');
    fclose(f)

    BIC = (-2 * LLtrace(length(LLtrace))) + (numPar * log(numDataPts));
    fname2 = [curDir '/BIC_' num2str(numClusters) '.mat'];
    f2 = fopen(fname2, 'at');
    cold2 = [num2str(numClusters) '\t' num2str(BIC, '%.6f')];
    fprintf(f2, cold2);
    fprintf(f2, '\n');
    fclose(f2)


    %%%%%%%%%%%%%%clustering script

    %%%start inference%%%
    disp(['Starting Inference...']);
    logl=0;
    marginal = zeros(numClusters, numDataPts);
    for (nd = 1:numDataPts)
        evidence = cell(numNodes, 1);
        evidence(2:numNodes,1) = data(2:numNodes,nd);
        engine = jtree_inf_engine(tan);
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
    outFile = [curDir '/classes-' num2str(numClusters) '.txt'];
    classe = [classes, marginal];
    display(size(classe));
    dlmwrite(outFile, classe, '\t');

end


exit

