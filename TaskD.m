%% Section D Clustering & Classification (note normalization as pre-processing)

%% 
% THIS CODE HAS BEEN MODIFIED FOR REPORT SUBMISSION. PLOTS HAVE BEEN 
% COMBINED TWO FIGURES (FOR PVT & ELECTRODE DATA). DECISION TREES WERE 
% GENERATE USING THE VIEW() COMMAND WHICH CANNOT BE ADDED TO FIGURES, SO
% THEY ALSO POP-UP.... SORRY :( ALSO NOTE VAIRABLES ARE CLEARED AT LINE 200 
% TO SHIFT FROM CLUSTERING ANALYSIS TO BAGGING ANALYSIS
%%

clear all;
close all;

load Data/F1_PVT.mat;

PVT = normalize( PVT , 1 );     %normalize data to stdev 1 and mean 0
numClusters = 6;

%% Hierarchical Clustering

% distance measure 
distEu = pdist(PVT , 'euclidean' );
distCB = pdist(PVT , 'cityblock' );

%links points according to closest distance from cluster average
ZEu = linkage(distEu , 'average' );
ZCB = linkage(distCB , 'average' );

%assigns each data point to a cluster number between 1 & MaxClust
c(:,1) = cluster( ZEu , 'Maxclust' , numClusters );
c(:,2) = cluster( ZCB , 'Maxclust' , numClusters );

%display clustering results
figure;
for j = 1 : 2
    
    subplot(3,2,j)
    for i = 1 : height(PVT)

        class = c(i,j);
        %choose plot color
        if class == 1
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'r' );
        elseif class == 2
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'g' );
        elseif class == 3
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'b' );
        elseif class == 4
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'm' );
        elseif class == 5
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'c'  );
        else
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'y'  );
        end 
        hold on;
    end
    ylabel( 'Vibrations (AU)' ) 
    xlabel( 'Pressure (AU)' ) 
    zlabel( 'Temperature (AU)' )
    if j == 1
        title( 'Hierarchical using Euclidean Distance' )
    else
        title( 'Hierarchical using City Block Distance' )
    end
end


%% K Means Clustering

%apply k means clustering using euclidean & city block distances
[idxEu,centEu] = kmeans(PVT,numClusters,'Distance','sqeuclidean');
[idxCB,centCB] = kmeans(PVT,numClusters,'Distance','cityblock');


for j = 3:4
    
    subplot(3,2,j);
    
    if j == 1 
        idx = idxEu;
    else
        idx = idxCB;
    end
    
    for i = 1 : height(PVT)

        if idx(i) == 1
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'r' );
        elseif idx(i) == 2
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'g' );
        elseif idx(i) == 3
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'b' );
        elseif idx(i) == 4
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'm' );
        elseif idx(i) == 5
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'c'  );
        else
            scatter3( PVT(i,1) , PVT(i,2) , PVT(i,3) , 'filled' , 'y'  );
        end 
        hold on;
    end

    for i = 1 : numClusters
        if j == 1 
            scatter3(centEu(i,1),centEu(i,2),centEu(i,3),'kx');
        else
            scatter3(centCB(i,1),centCB(i,2),centCB(i,3),'kx');
        end
        hold on;
    end
    
    if j == 3
        title( 'k-means using Euclidean Distance' )
    else
        title( 'k-means using City Block Distance' )
    end

    ylabel( 'Vibrations (AU)' ) 
    xlabel( 'Pressure (AU)' ) 
    zlabel( 'Temperature (AU)' )
end

%% Principal Component Analysis

featureAvg = mean( PVT , 1 );                       %finds average of each feature
averageMatrix = ones( size(PVT) ).* featureAvg;     %creates average matrix
standPVT = PVT - featureAvg;                        %standardizes data matrix
covarinaceMatrix = cov(standPVT);                   %covariance matrix calculation

[ U S V ] = svd( covarinaceMatrix , 'econ' );       %find eigenvalues and eigenvectors

% PVT reproject
for i = 1 : height(standPVT)
   
    x = V(:,1)' * standPVT(i,:)';
    y = V(:,2)' * standPVT(i,:)';
    z = V(:,3)' * standPVT(i,:)';
    
    rPVT(i,:) = [ x , y , z ];  
end

rPVT = normalize(rPVT);

%% K Means Clustering Reprojected Data

[idxEu,centEu] = kmeans(rPVT,numClusters,'Distance','sqeuclidean');
[idxCB,centCB] = kmeans(rPVT,numClusters,'Distance','cityblock');


for j = 5:6
    
    subplot(3,2,j);
    
    if j == 1 
        idx = idxEu;
    else
        idx = idxCB;
    end
    
    for i = 1 : height(PVT)

        if idx(i) == 1
            scatter3( rPVT(i,1) , rPVT(i,2) , rPVT(i,3) , 'filled' , 'r' );
        elseif idx(i) == 2
            scatter3( rPVT(i,1) , rPVT(i,2) , rPVT(i,3) , 'filled' , 'g' );
        elseif idx(i) == 3
            scatter3( rPVT(i,1) , rPVT(i,2) , rPVT(i,3) , 'filled' , 'b' );
        elseif idx(i) == 4
            scatter3( rPVT(i,1) , rPVT(i,2) , rPVT(i,3) , 'filled' , 'm' );
        elseif idx(i) == 5
            scatter3( rPVT(i,1) , rPVT(i,2) , rPVT(i,3) , 'filled' , 'c'  );
        else
            scatter3( rPVT(i,1) , rPVT(i,2) , rPVT(i,3) , 'filled' , 'y'  );
        end 
        hold on;
    end

    for i = 1 : numClusters
        if j == 1 
            scatter3(centEu(i,1),centEu(i,2),centEu(i,3),'kx');
        else
            scatter3(centCB(i,1),centCB(i,2),centCB(i,3),'kx');
        end
        hold on;
    end
    
    if j == 5
        title( 'PCA & k-means using Euclidean Distance' )
    else
        title( 'PCA & k-means using City Block Distance' )
    end

    ylabel( 'Vibrations (AU)' ) 
    xlabel( 'Pressure (AU)' ) 
    zlabel( 'Temperature (AU)' )
end

%% Bagging

clear all;
load Data/F1_Elec;

%% Load data & applpy PCA

featureAvg = mean( Electrode , 1 );                     %finds average of each feature
averageMatrix = ones( size(Electrode) ).* featureAvg;   %creates average matrix
standElec = Electrode - averageMatrix;                  %standardizes data matrix
covMatElec = cov(standElec);                            %find covariance matrix

[ U S V ] = svd(covMatElec , 'econ' );                  %find eigenvalues and eigenvectors

numTrials = 10;

for i = 1 : numTrials : length(standElec)
    %reproject data using 3 largest PCs
    rPVT(i:i+numTrials-1,1) = ( V(:,1)' * standElec(i:i+numTrials-1,:)' )' ;
    rPVT(i:i+numTrials-1,2) = ( V(:,2)' * standElec(i:i+numTrials-1,:)' )' ;
    rPVT(i:i+numTrials-1,3) = ( V(:,3)' * standElec(i:i+numTrials-1,:)' )' ;
end

%adds label to data
class = [ ones(numTrials,1) ; 2*ones(numTrials,1) ; 3*ones(numTrials,1) ; 4*ones(numTrials,1) ; 5*ones(numTrials,1) ; 6*ones(numTrials,1) ];
rPVT = cat(2,rPVT,class);

%% Split data

trainingSplit = 0.6 * height(rPVT);
testingSplit = 0.4 * height(rPVT);

if trainingSplit+testingSplit ~= height(rPVT)
    error("Incorrect Training & Testing Splits")
end

ix = randperm(height(rPVT));    %randomly pick trials

%place trails in testing or training
trainPVT = rPVT( ix(1:trainingSplit) , : );
testPVT = rPVT( ix(1:testingSplit) , : );

%% Bagging

B = TreeBagger( 8 , trainPVT(:,1:3) , trainPVT(:,4) , 'Method' , 'Classification' , 'OOBPrediction', 'On' );

%visualize to choose number of trees (8-10)
figure;
subplot(1,2,1)
errorOOB = oobError(B);
plot(errorOOB)
xlabel ('Number of grown trees' , 'FontSize' , 16 );
ylabel ( 'Out-of-bag classification error' , 'FontSize' , 16 );

%displays some trees
view(B.Trees{1},'Mode','graph')
%view(B.Trees{5},'Mode','graph') %commented out for report to reduce number of windows that open up

%get class predictions for test data using trained bagging algorithm
pred = predict( B , testPVT(:,1:3) );

%convert predictions to vector
pred = cell2mat(pred);
prediction = zeros(length(pred),1);
for i = 1 : length(pred)
    prediction(i) = str2num(pred(i));
end

subplot(1,2,2)
cm = confusionchart( testPVT(:,4) , prediction );




