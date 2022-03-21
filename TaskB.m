%% Section B: Principal Component Analysis

%% 
% THIS CODE HAS BEEN MODIFIED FOR REPORT SUBMISSION. PLOTS HAVE BEEN 
% COMBINED ONE FIGURE. AS IT IS VERY FULL LEGENDS HAVE BEEN COMMENTED OUT.
% ALSO NOTE VAIRABLES ARE CLEARED AT LINE 104 TO SHIFT FROM PVT ANALYSIS
% TO ELECTRODE ANALYSIS
%%

close all; clear all;

load Data/F1_PVT

%% 1a Setup

featureAvg = mean( PVT , 1 );                       %finds average of each feature
averageMatrix = ones( size(PVT) ).* featureAvg;     %creates average matrix

standPVT = PVT - averageMatrix;                     %standardizes data matrix

covarinaceMatrix = cov(standPVT);                   %finds the covariace matrix

[ U S V ] = svd( covarinaceMatrix , 'econ' );       %find eigenvalues and eigenvectors

%% 1b Plot

numTrials = 10; %number of trials per object

figure;
subplot(2,4,1)
for i = 1 : numTrials : height(PVT)
        %plot standardized data
        scatter3( standPVT(i:i+numTrials-1,1) , standPVT(i:i+numTrials-1,2) , standPVT(i:i+numTrials-1,3) , 'filled' );
        hold on;
end

%plot principal components
%PC1 (note s1 is determined empirically it scales the eigenvector for visulization purposes)
s1 = -200 : 200;
plot3( V(1,1)*s1 , V(2,1)*s1 , V(3,1)*s1 , 'LineWidth' , 1.5 );
hold on;
%PC 2
s2 = -50 : 50;
plot3( V(1,2)*s2 , V(2,2)*s2 , V(3,2)*s2 , 'LineWidth' , 1.5 );
%PC 2
s3 = -20 : 20;
plot3( V(1,3)*s3 , V(2,3)*s3 , V(3,3)*s3 , 'LineWidth' , 1.5 );
%commented out for report submission
%legend( 'Acrylic' , 'Black Foam' , 'Car Sponge' , 'Flour Sack' , 'Kitchen Sponge' , 'Steel Vase' , 'Principal Component 1' , 'Principal Component 2' , 'Principal Component 3' , 'FontSize' , 16 );
xlabel( 'Standarsied Vibrations (AU)' , 'FontSize' , 16 ) 
ylabel( 'Standarsied Pressure (AU)' , 'FontSize' , 16 ) 
zlabel( 'Standarsied Temperature (AU)' , 'FontSize' , 16 )

%% 1c Reduce data to 2D

%replot data along 2 largest PC axes
subplot(2,4,2)
for i = 1 : numTrials : length(standPVT)
    
    %reprojet along principal axes
    PC1 = V(:,1)' * standPVT(i:i+numTrials-1,:)' ;
    PC2 = V(:,2)' * standPVT(i:i+numTrials-1,:)' ;
    
    scatter( PC1 , PC2 , 'filled' );
    hold on;
end
%commented out for report submission
%legend( 'Acrylic' , 'Black Foam' , 'Car Sponge' , 'Flour Sack' , 'Kitchen Sponge' , 'Steel Vase' , 'FontSize' , 16 );
xlabel( 'Principal Component 1' , 'FontSize' , 16 )
ylabel( 'Principal Component 2' , 'FontSize' , 16 )

%% 1d 1D lines

y = zeros(numTrials,1);

for i = 1 : numTrials : length(standPVT)
    
    %find data points along PC axes
    PC1 = V(:,1)' * standPVT(i:i+numTrials-1,:)' ;
    PC2 = V(:,2)' * standPVT(i:i+numTrials-1,:)' ;
    PC3 = V(:,3)' * standPVT(i:i+numTrials-1,:)' ;

    subplot(2,4,3);
    plot(PC1 , y , '.' , 'MarkerSize' , 10);
    title( 'Principal Component 1' , 'FontSize' , 16 );
    hold on;
    
    subplot(2,4,4);
    plot(PC2 , y , '.' , 'MarkerSize' , 10);
    title( 'Principal Component 2' , 'FontSize' , 16 );
    hold on;

    subplot(2,4,5);
    plot(PC3 , y , '.' , 'MarkerSize' , 10);
    title( 'Principal Component 3' , 'FontSize' , 16);
    hold on;
    
end
%commented out for report submission
%legend( 'Acrylic' , 'Black Foam' , 'Car Sponge' , 'Flour Sack' , 'Kitchen Sponge' , 'Steel Vase' , 'FontSize' , 16 );

%% 2a Scree Plot for electrodes

clear all;

load Data/F1_Elec

featureAvg = mean( Electrode , 1 );                     %finds average of each feature
averageMatrix = ones( size(Electrode) ).* featureAvg;   %creates average matrix

standElec = Electrode - averageMatrix;                  %standardizes data matrix

covMatElec = cov(standElec);

[ U S V ] = svd(covMatElec , 'econ' );                   %find eigenvalues and eigenvectors

eigenvalues = max(S);   %creates a vector of eigenvalues

%scree plot
subplot(2,4,6)
plot ( [eigenvalues] , 'b--o' )
xlabel( 'Principal Components' , 'FontSize' , 16 )
ylabel( 'Variance (or eigenvalue)' , 'FontSize' , 16 )


% 2b 3 Largest PCs Visualization

numTrials = 10;

subplot(2,4,7)
for i = 1 : numTrials : length(standElec)
    
    PC1 = V(:,1)' * standElec(i:i+numTrials-1,:)' ;
    PC2 = V(:,2)' * standElec(i:i+numTrials-1,:)' ;
    PC3 = V(:,3)' * standElec(i:i+numTrials-1,:)' ;
    
    scatter3( PC1 , PC2 , PC3 , 'filled' );
    hold on;
end
%commented out for report submission
%legend( 'Acrylic' , 'Black Foam' , 'Car Sponge' , 'Flour Sack' , 'Kitchen Sponge' , 'Steel Vase' , 'FontSize' , 12 );
xlabel( 'Principal Component 1' , 'FontSize' , 16)
ylabel( 'Principal Component 2' , 'FontSize' , 16)
zlabel( 'Principal Component 3' , 'FontSize' , 16)

%% Choosing electrodes

%pick most useful electrodes
for i = 1:19
    inter(:,i) = ((U(:,i)).^2).*S(i,i);
end
for i = 1:19
    sumv(i) = sum(inter(:,i));
end

subplot(2,4,8)
bar( [log10(sumv)] )
xlabel( 'Electrode Number' , 'Fontsize' , 16 )
ylabel( 'log10( Measure of data captured) ' , 'FontSize' , 16 )

