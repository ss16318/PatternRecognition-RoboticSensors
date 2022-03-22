%% Task C

%% 
% THIS CODE HAS BEEN MODIFIED FOR REPORT SUBMISSION. PLOTS HAVE BEEN 
% COMBINED INTO ONE FIGURE. AS IT IS VERY FULL, LEGENDS HAVE BEEN COMMENTED 
% OUT.
%%

clear all; close all;

load Data/F1_PVT

%% LDA for PVT

numTrials = 10; %number of trials per object

%% Step 1 Standardize

%create array of objects for comparison
objects(1:numTrials,:) = PVT( 11:20 , :);
objects(numTrials+1:numTrials+10,:) = PVT( 21:30 , :);

%standardizes data to have zero mean and stdev of 1
normObjs = normalize(objects , 1);

%% Step 2 Find Data Means

avgObjs(1,:) = mean( normObjs(1:numTrials,:) ,1 );
avgObjs(2,:) = mean( normObjs(numTrials+1:numTrials*2,:) , 1 );

%% Step 3 Find Between class scatter matrix (check)
Sb = zeros(2,2,3);

%find differences between class averages
difMeans = ( avgObjs(1,:) - avgObjs(2,:) )';
%pressure vs vibration difference
DIF1 = [ difMeans(1) ; difMeans(2) ];
Sb(:,:,1) = DIF1 * DIF1' ; 
%pressure vs temperature difference
DIF2 = [ difMeans(1) ; difMeans(3) ];
Sb(:,:,2) = DIF2 * DIF2' ; 
%temperature vs vibration difference
DIF3 = [ difMeans(3) ; difMeans(2) ];
Sb(:,:,3) = DIF3 * DIF3' ; 

%% Find Within class scatter matrix (check)

for i = 1 : 3
    %split data for each object
    X1 = normObjs(1:numTrials,i);
    X2 = normObjs(numTrials+1:2*numTrials,i);
    %find class centroids
    avg1 = avgObjs(1,i) * ones( size(X1) );
    avg2 = avgObjs(2,i) * ones( size(X2) );
    %find difference between class data and class centroid
    dif(1,:,i) = (X1 - avg1)';
    dif(2,:,i) = (X2 - avg2)'; 
    
end

Sw = zeros(2,2,3);  %initialize within class scatter matrices

%pressure vs vibration
dif1 = [ dif(1,:,1) ; dif(1,:,2) ];
dif2 = [ dif(2,:,1) ; dif(2,:,2) ];
Sw(:,:,1) = (dif1 * dif1') + (dif2 * dif2') ;

%pressure vs temperature
dif1 = [ dif(1,:,1) ; dif(1,:,3) ];
dif2 = [ dif(2,:,1) ; dif(2,:,3) ];
Sw(:,:,2) = (dif1 * dif1') + (dif2 * dif2') ;

%temperature vs vibration
dif1 = [ dif(1,:,3) ; dif(1,:,2) ];
dif2 = [ dif(2,:,3) ; dif(2,:,2) ];
Sw(:,:,3) = (dif1 * dif1') + (dif2 * dif2') ;

%% Eigenmodes

for i = 1 : size(Sw,3)
   
    %create data matrix of scatter matrices
    A = inv(Sw(:,:,i)) * Sb(:,:,i);
    %find eigenmodes
    [V,D] = eig(A);
    %save eigenmodes
    eigV(:,:,i) = V;
    eigD(:,:,i) = D;
    
    [ Dmax(i) , I ] = max(max(D)); %find largest eigenvalue
    LDAV(:,i) = V(:,I);            %save corresponding eigenvector
    
end

%% Plots for a

%plot pressure vs vibration
figure;
subplot(2,5,1)
for i = 1 : numTrials : height(normObjs)
    %determine point colors    
    if i <= 10
        c = '#D95319';
    else
        c = '#EDB120';
    end
    
    scatter( normObjs(i:i+numTrials-1,2) , normObjs(i:i+numTrials-1,1) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end
% plots means of data groups
plot( avgObjs(1,2) , avgObjs(1,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#D95319' );
plot( avgObjs(2,2) , avgObjs(2,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#EDB120' );

% plot LDA lines
s = -2 : 2;
plot( LDAV(2,1)*s , LDAV(1,1)*s , 'LineWidth' , 1.5 );

xlabel( 'Vibration (AU)' )
ylabel( 'Pressure (AU)' )
title( 'Pressure against Vibration' ) %, 'FontSize' , 16 )

%plot pressure vs temperature
subplot(2,5,2)
for i = 1 : numTrials : height(normObjs)
        
    if i <= 10
        c = '#D95319';
    else
        c = '#EDB120';
    end
    
    scatter( normObjs(i:i+numTrials-1,3) , normObjs(i:i+numTrials-1,1) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end
% plots means of data groups
plot( avgObjs(1,3) , avgObjs(1,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#D95319' );
plot( avgObjs(2,3) , avgObjs(2,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#EDB120');

% Plot LDA lines
plot( LDAV(2,2)*s , LDAV(1,2)*s , 'LineWidth' , 1.5 );

xlabel( 'Temperature (AU)' )
ylabel( 'Pressure (AU)' )
title( 'Pressure against Temperature' ) %, 'FontSize' , 16 )

%plot temperature vs vibration
subplot(2,5,3)
for i = 1 : numTrials : height(normObjs)
        
    if i <= 10
        c = '#D95319';
    else
        c = '#EDB120';
    end
    
    scatter( normObjs(i:i+numTrials-1,2) , normObjs(i:i+numTrials-1,3) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end

% plots means of data groups
plot( avgObjs(1,2) , avgObjs(1,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#D95319' );
plot( avgObjs(2,2) , avgObjs(2,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#EDB120' );

% Plot LDA lines
plot( LDAV(2,3)*s , LDAV(1,3)*s , 'LineWidth' , 1.5 );

ylabel( 'Temperature (AU)')
xlabel( 'Vibration (AU)' )
title( 'Temperature against Vibration' ) %, 'FontSize' , 16 )
%the line below was commented out for the code submission
%legend( 'Black Foam' , 'Car Sponge' , 'Mean of Black Foam' , 'Mean of Car Sponge' , 'LDA Line' , 'FontSize' , 16 );

%% LDA 3D reprojection

% Between class scatter matrix
Sb3 = difMeans * difMeans';

% Within class scatter matrix
DIF1 = [ dif(1,:,1) ; dif(1,:,2) ; dif(1,:,3) ];
DIF2 = [ dif(2,:,1) ; dif(2,:,2) ; dif(2,:,3) ];

Sw3 = ( DIF1 * DIF1' ) + ( DIF2 * DIF2' );

%create data matrix of scatter matrices
A3 = inv(Sw3) * Sb3;
%find eigenmodes
[V3,D3] = eig(A3);

%finds eigenmode w/ biggest eigenvalue
[ maxD3 I ] = max(max(D3));
V3LDA = V3(:,I);

%find eigenmode w/ second biggest eigenvalue
D3(I,I) = 0;    %sets largest eigenvalue to 0
[ D3next I2 ] = max(max(D3));
V3LDA2 = V3(:,I2);

%% plot in 3D

%plot PVT data
subplot(2,5,4)
for i = 1 : numTrials : height(normObjs)
        
    if i <= 10
        c = '#D95319';
    else
        c = '#EDB120';
    end
    
    scatter3( normObjs(i:i+numTrials-1,1) , normObjs(i:i+numTrials-1,2) , normObjs(i:i+numTrials-1,3) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end
% plot means
plot3( avgObjs(1,1) , avgObjs(1,2) , avgObjs(1,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#D95319' );
plot3( avgObjs(2,1) , avgObjs(2,2) , avgObjs(2,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#EDB120' );

%plot LDA lines
plot3( V3LDA(1)*s , V3LDA(2)*s , V3LDA(3)*s , 'LineWidth' , 1.5 );

% the line below was commented out for the code submission
%legend( 'Black Foam' , 'Car Sponge' , 'Mean of Black Foam' , 'Mean of Car Sponge' , 'LDA 1' );
xlabel( 'Vibrations (AU)' ) 
ylabel( 'Pressure (AU)' ) 
zlabel( 'Temperature (AU)' )

%% Plot LDA reduction

y = zeros( length(normObjs)/2 ,1 );

for i = 1 : numTrials : length(normObjs)
        
    %reprojection
    X1(:) = V3LDA' * normObjs(i:i+numTrials-1,:)';
    X2(:) = V3LDA2' * normObjs(i:i+numTrials-1,:)';
    
    if i <= 10
        c = '#D95319';
    else
        c = '#EDB120';
    end
    
    subplot(2,5,5)
    plot(X1 , y , '.' , 'MarkerSize' , 30 , 'Color' , c );
    title( 'LDA Reduction to 1D' ) %, 'FontSize' , 16 )
    hold on;
end


%% Acrylic & Steel Vase


%% LDA for PVT

%% Step 1 Standardize:

%create array of objects for comparison
objects(1:numTrials,:) = PVT( 1:10 , :);
objects(numTrials+1:numTrials+10,:) = PVT( 51:60 , :);

%standardizes data to have zero mean and stdev of 1
normObjs = normalize(objects , 1);

%% Step 2 Find Data Means

avgObjs(1,:) = mean( normObjs(1:numTrials,:) ,1 );
avgObjs(2,:) = mean( normObjs(numTrials+1:numTrials*2,:) , 1 );

%% Step 3 Find Between class scatter matrix (check)

Sb = zeros(2,2,3);

%find differences between class averages
difMeans = ( avgObjs(1,:) - avgObjs(2,:) )';
%pressure vs vibration difference
DIF1 = [ difMeans(1) ; difMeans(2) ];
Sb(:,:,1) = DIF1 * DIF1' ; 
%pressure vs temperature difference
DIF2 = [ difMeans(1) ; difMeans(3) ];
Sb(:,:,2) = DIF2 * DIF2' ; 
%temperature vs vibration difference
DIF3 = [ difMeans(3) ; difMeans(2) ];
Sb(:,:,3) = DIF3 * DIF3' ; 

%% Find Within class scatter matrix (check)

for i = 1 : 3
    %splits object data
    X1 = normObjs(1:numTrials,i);
    X2 = normObjs(numTrials+1:2*numTrials,i);
    %gets object centroids
    avg1 = avgObjs(1,i) * ones( size(X1) );
    avg2 = avgObjs(2,i) * ones( size(X2) );
    
    dif(1,:,i) = (X1 - avg1)';
    dif(2,:,i) = (X2 - avg2)'; 
    
end

Sw = zeros(2,2,3);

%pressure vs vibration
dif1 = [ dif(1,:,1) ; dif(1,:,2) ];
dif2 = [ dif(2,:,1) ; dif(2,:,2) ];
Sw(:,:,1) = (dif1 * dif1') + (dif2 * dif2') ;

%pressure vs temperature
dif1 = [ dif(1,:,1) ; dif(1,:,3) ];
dif2 = [ dif(2,:,1) ; dif(2,:,3) ];
Sw(:,:,2) = (dif1 * dif1') + (dif2 * dif2') ;

%temperature vs vibration
dif1 = [ dif(1,:,3) ; dif(1,:,2) ];
dif2 = [ dif(2,:,3) ; dif(2,:,2) ];
Sw(:,:,3) = (dif1 * dif1') + (dif2 * dif2') ;

%% Eigenmodes

for i = 1 : size(Sw,3)
   
    %create data matrix of scatter matrices
    A = inv(Sw(:,:,i)) * Sb(:,:,i);
    %find eigenmodes
    [V,D] = eig(A);
    %saves eigenmodes
    eigV(:,:,i) = V;
    eigD(:,:,i) = D;
    
    [ Dmax(i) , I ] = max(max(D));  %finds largest eigenvalue
    LDAV(:,i) = V(:,I);             %saves corresponding eigenvalue
    
end

%% Plots for a

%plot pressure vs vibration
subplot(2,5,6)
for i = 1 : numTrials : height(normObjs)
        
    if i <= 10
        c = '#0000FF';
    else
        c = '#00FFFF';
    end
    
    scatter( normObjs(i:i+numTrials-1,2) , normObjs(i:i+numTrials-1,1) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end
% plots means of data groups
plot( avgObjs(1,2) , avgObjs(1,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#0000FF' );
plot( avgObjs(2,2) , avgObjs(2,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#00FFFF' );

% plot LDA lines
s = -2 : 2;
plot( LDAV(2,1)*s , LDAV(1,1)*s , 'LineWidth' , 1.5 );

% the next line was commented out for the report
%legend( 'Acrylic' , 'Steel Vase' , 'Mean of Acrylic' , 'Mean of Steel Vase' , 'LDA Line' );
xlabel( 'Vibration (AU)' )
ylabel( 'Pressure (AU)'  )
title( 'Pressure against Vibration')

%plot pressure vs temperature
subplot(2,5,7)
for i = 1 : numTrials : height(normObjs)
        
    if i <= 10
        c = '#0000FF';
    else
        c = '#00FFFF';
    end
    
    scatter( normObjs(i:i+numTrials-1,3) , normObjs(i:i+numTrials-1,1) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end
% plots means of data groups
plot( avgObjs(1,3) , avgObjs(1,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#0000FF' );
plot( avgObjs(2,3) , avgObjs(2,1) , '*' , 'MarkerSize' , 10 , 'Color' , '#00FFFF');

% Plot LDA lines
plot( LDAV(2,2)*s , LDAV(1,2)*s , 'LineWidth' , 1.5 );

xlabel( 'Temperature (AU)')
ylabel( 'Pressure (AU)' )
title( 'Pressure against Temperature')

%plot temperature vs vibration
subplot(2,5,8)
for i = 1 : numTrials : height(normObjs)
        
    if i <= 10
        c = '#0000FF';
    else
        c = '#00FFFF';
    end
    
    scatter( normObjs(i:i+numTrials-1,2) , normObjs(i:i+numTrials-1,3) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end

% plots means of data groups
plot( avgObjs(1,2) , avgObjs(1,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#0000FF' );
plot( avgObjs(2,2) , avgObjs(2,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#00FFFF' );

% Plot LDA lines
plot( LDAV(2,3)*s , LDAV(1,3)*s , 'LineWidth' , 1.5 );

ylabel( 'Temperature (AU)')
xlabel( 'Vibration (AU)' )
title( 'Temperature against Vibration')

%% LDA 3D reprojection

% Between class scatter matrix
Sb3 = difMeans * difMeans';

% Within class scatter matrix

DIF1 = [ dif(1,:,1) ; dif(1,:,2) ; dif(1,:,3) ];
DIF2 = [ dif(2,:,1) ; dif(2,:,2) ; dif(2,:,3) ];

Sw3 = ( DIF1 * DIF1' ) + ( DIF2 * DIF2' );

%create data matrix of scatter matrices
A3 = inv(Sw3) * Sb3;
%find eigenmodes (issues starts)
[V3,D3] = eig(A3);

%finds eigenmode w/ biggest eigenvalue
[ maxD3 , I ] = max(max(D3));
V3LDA = V3(:,I);

%find eigenmode w/ second biggest eigenvalue
D3(I,I) = 0;    %sets largest eigenvalue to 0
[ D3next I2 ] = max(max(D3));
V3LDA2 = V3(:,I2);

%% plot in 3D

%plot PVT data
subplot(2,5,9)
for i = 1 : numTrials : height(normObjs)
        
    if i <= 10
        c = '#0000FF';
    else
        c = '#00FFFF';
    end
    
    scatter3( normObjs(i:i+numTrials-1,1) , normObjs(i:i+numTrials-1,2) , normObjs(i:i+numTrials-1,3) , 'filled' , 'MarkerFaceColor' , c );
    hold on;
end
% plot means
plot3( avgObjs(1,1) , avgObjs(1,2) , avgObjs(1,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#0000FF' );
plot3( avgObjs(2,1) , avgObjs(2,2) , avgObjs(2,3) , '*' , 'MarkerSize' , 10 , 'Color' , '#00FFFF' );

%plot LDA lines
plot3( V3LDA(1)*s , V3LDA(2)*s , V3LDA(3)*s , 'LineWidth' , 1.5 );
plot3( V3LDA2(1)*s , V3LDA2(2)*s , V3LDA2(3)*s , 'LineWidth' , 1.5 );

%the line below was commented out for the code submission
%legend( 'Acrylic' , 'Steel Vase' , 'Mean of Acrylic' , 'Mean of Steel Vase' , 'LDA 1' , 'LDA 2' );
xlabel( 'Vibrations (AU)' ) 
ylabel( 'Pressure (AU)' ) 
zlabel( 'Temperature (AU)' )

%% Plot LDA reduction

y = zeros( length(normObjs)/2 ,1 );
subplot(2,5,10);
for i = 1 : numTrials : length(normObjs)
        
    %reprojection
    X1(:) = V3LDA' * normObjs(i:i+numTrials-1,:)';
    X2(:) = V3LDA2' * normObjs(i:i+numTrials-1,:)';
    
    if i <= 10
        c = '#0000FF';
    else
        c = '#00FFFF';
    end
    
    plot(X1 , y , '.' , 'MarkerSize' , 30 , 'Color' , c );
    title( 'LDA Reduction to 1D' )
    hold on;
end
