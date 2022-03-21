%% Task A

%% 
% THIS CODE HAS BEEN MODIFIED FOR REPORT SUBMISSION. NOTE VAIRABLES ARE 
% CLEARED AT LINE 44 TO SHIFT FROM RAW DATA ANALYSIS TO DATA PREPARATION
% AND AT LINE 315 TO SHIFT FROM DATA PREP TO TIMESTEP VISUALIZATION
%%

clear all;
close all;

%% Visualising raw data

files = dir('*.mat');
load( files(45).name );     %objects 5 & 45 were shown in the report

time = [ 1 : 1000 ];

figure;
subplot(2,2,1)
for i = 1 : height(F1Electrodes)
    plot( time , F1Electrodes(i,:) );
    title('Electrodes')
    hold on;
end

subplot(2,2,2)
plot( time , F1pac(2,:) );
title('Vibrations');
ylim( [1800 2100] );

subplot( 2 ,2 ,3 );
plot( time , F1pdc );
title('Pressure');
ylim( [1050 1320] );

subplot( 2 ,2 ,4 );
plot( time , F1tdc );
title('Temperature');
ylim( [1920 1950] );

%% Choosing specific timestep and finer algorithm

clear all;

candidatestep = zeros(12,1);

%preallocate memory
F0Electrodes3D = zeros(60,19,1000);
F0pac3D = zeros(60,1,1000);
F0pdc3D = zeros(60,1,1000);
F0tdc3D = zeros(60,1,1000);
F1Electrodes3D = zeros(60,19,1000);
F1pac3D = zeros(60,1,1000);
F1pdc3D = zeros(60,1,1000);
F1tdc3D = zeros(60,1,1000);

%averaging in time for each trial
F0Electrodesavg = zeros(10,19);
F0pacavg = zeros(10,1);
F0pdcavg = zeros(10,1);
F0tdcavg = zeros(10,1);
F1Electrodesavg = zeros(10,19);
F1pacavg = zeros(10,1);
F1pdcavg = zeros(10,1);
F1tdcavg = zeros(10,1);

%load raw data files
files = dir('*.mat');

T = 100;    %maxTime

for k = 1 : size(files,1)
    
    load( files(k).name );  %load each file
    
    for i = 1 : T    %loop through each time instant
        
        F0Electrodes3D(k,:,i) = F0Electrodes(:,i);
        F0pac3D(k,:,i) = F0pac(2,i);
        F0pdc3D(k,:,i) = F0pdc(:,i);
        F0tdc3D(k,:,i) = F0tdc(:,i);
        F1Electrodes3D(k,:,i) = F1Electrodes(:,i);
        F1pac3D(k,:,i) = F1pac(2,i);
        F1pdc3D(k,:,i) = F1pdc(:,i);
        F1tdc3D(k,:,i) = F1tdc(:,i);  
    end
end

%store standard deviations of objects at each time instant
F0Electrodesstd = zeros(60,T,19);
F1Electrodesstd = zeros(60,T,19);
F0pacstd = zeros(60,T);
F0pdcstd = zeros(60,T);
F0tdcstd = zeros(60,T);
F1pacstd = zeros(60,T);
F1pdcstd = zeros(60,T);
F1tdcstd = zeros(60,T);

idx = [1 10 11 20 21 30 31 40 41 50 51 60]; %take indexes of change of object

for u = 1:60
   
    %identify object
    if u < 11
        n = 1;
    end
    if (u > 10) && (u < 21)
        n = 2;
    end
    if (u > 20) && (u < 31)
        n = 3;
    end
    if (u > 30) && (u < 41)
        n = 4;
    end
    if (u > 40) && (u < 51)
        n = 5;
    end
    if u > 50
        n = 6;
    end
    
    for i = 1 : T       %loops through each time instant & finds std across object class
        
        for j = 1:19
            F0Electrodesstd(u,i,j) = std( F0Electrodes3D( idx(2*n-1):idx(2*n),j,i ) );
            F1Electrodesstd(u,i,j) = std( F1Electrodes3D( idx(2*n-1):idx(2*n),j,i ) );
        end
        
        F0pacstd(u,i,1) = std( F0pac3D( idx(2*n-1):idx(2*n),i ) );
        F0pdcstd(u,i,1) = std( F0pdc3D( idx(2*n-1):idx(2*n),i ) );
        F0tdcstd(u,i,1) = std( F0tdc3D( idx(2*n-1):idx(2*n),i ) );
        F1pacstd(u,i,1) = std( F1pac3D( idx(2*n-1):idx(2*n),i ) );
        F1pdcstd(u,i,1) = std( F1pdc3D( idx(2*n-1):idx(2*n),i ) );
        F1tdcstd(u,i,1) = std( F1tdc3D( idx(2*n-1):idx(2*n),i ) );
    end
end

%find max value of 8 attributes for normalization
for n = 1 : 6
    for i = 1 : T
        for c = 1 : 10    %loops through the 10 trials for each object
            
            v1(c) = max( F0Electrodesstd((n-1)*10+c,i,:) );   %max of 19 sensors
            v2(c) = max( F1Electrodesstd((n-1)*10+c,i,:) );
            v3(c) = F0pacstd( (n-1)*10+c,i );
            v4(c) = F0pdcstd( (n-1)*10+c,i );
            v5(c) = F0tdcstd( (n-1)*10+c,i );
            v6(c) = F1pacstd( (n-1)*10+c,i );
            v7(c) = F1pdcstd( (n-1)*10+c,i );
            v8(c) = F1tdcstd( (n-1)*10+c,i );
        end
        
        F0Electrodesstdmax(n,i) = max(v1);
        F1Electrodesstdmax(n,i) = max(v2);
        F0pacstdmax(n,i) = max(v3);
        F0pdcstdmax(n,i) = max(v4);
        F0tdcstdmax(n,i) = max(v5);
        F1pacstdmax(n,i) = max(v6);
        F1pdcstdmax(n,i) = max(v7);
        F1tdcstdmax(n,i) = max(v8);
    end
end

for n = 1:6
    E0(n) = max(F0Electrodesstdmax(n,:));
    pac0(n) = max(F0pacstdmax(n,:));
    pdc0(n) = max(F0pdcstdmax(n,:));
    tdc0(n) = max(F0tdcstdmax(n,:));
    E1(n) = max(F1Electrodesstdmax(n,:));
    pac1(n) = max(F1pacstdmax(n,:));
    pdc1(n) = max(F1pdcstdmax(n,:));
    tdc1(n) = max(F1tdcstdmax(n,:));
end

stddev1 = zeros(6,T);
%find timesteps at which 10 trials are most similar in each object for F1
for u = 1:60
    
    if u < 11
        n = 1;
    end
    if (u > 10) && (u < 21)
        n = 2;
    end
    if (u > 20) && (u < 31)
        n = 3;
    end
    if (u > 30) && (u < 41)
        n = 4;
    end
    if (u > 40) && (u < 51)
        n = 5;
    end
    if u > 50
        n = 6;
    end
    
    for i = 1:T
        esection1 = 0;
        for j = 1:19
            esection1 = esection1 + F1Electrodesstd(u,i,j);
        end
        esection = esection1/(19*E1(n));
        asection = F1pacstd(u,i)/pac1(n);
        dsection = F1pdcstd(u,i)/pdc1(n);
        tsection = F1tdcstd(u,i)/tdc1(n);
        stddev1(n,i) = esection + asection + dsection + tsection;
    end
end

stddev0 = zeros(6,T);
%find timesteps at which 10 trials are most similar in each object for F0
for u = 1:60
    
    if u < 11
        n = 1;
    end
    if (u > 10) && (u < 21)
        n = 2;
    end
    if (u > 20) && (u < 31)
        n = 3;
    end
    if (u > 30) && (u < 41)
        n = 4;
    end
    if (u > 40) && (u < 51)
        n = 5;
    end
    if u > 50
        n = 6;
    end
    
    for i = 1:T
        esection0 = 0;
        for j = 1:19
            esection0 = esection0 + F0Electrodesstd(u,i,j);
        end
        esection = esection0/19/E0(n);
        asection = F0pacstd(u,i,1)/pac0(n);
        dsection = F0pdcstd(u,i,1)/pdc0(n);
        tsection = F0tdcstd(u,i,1)/tdc0(n);
        stddev0(n,i) = esection + asection + dsection + tsection;
    end
end

x0 = zeros(1,6);
for n = 1:6
    minstddev = min(stddev0(n,:));
    [x0(n),candidatestep(n,1)] = find(stddev0(n,:)==minstddev);
end

x1 = zeros(1,6);
for n = 1:6
    minstddev = min(stddev1(n,:));
    [x1(n),candidatestep(n+6,1)] = find(stddev1(n,:)==minstddev);
end

%find the 1 timestep in candidatestep where object separation is largest
candidatstddev = zeros(12,1);
for n = 1:6
    candidatstddev(n,1) = std(stddev0(:,candidatestep(n)));
end
for n = 7:12
    candidatstddev(n,1) = std(stddev1(:,candidatestep(n)));
end
%find maximum interclass stdev
[~, I] = max(candidatstddev);

Time = candidatestep(I);    %finds best timestep

%determines finger to use
if I <= 6
    finger = 0;
else
    finger = 1;
end

%% Creating .mat file

%remove subfolder called Data and its files as new data will be saved here
if isfolder( 'Data')
   rmdir('Data', 's')
end

files = dir('*.mat');   %load file names

for i = 1 : length(files)
    
    load( files(i).name ); %load each trial
    
    if finger == 0
        %store PVT data for each trial (note cols are pressure, vibration & temperature)
        PVT(i,:) = [ F0pdc(Time)  F0pac(2,Time) F0tdc(Time) ];
        %store electrde data for each trial
        Electrode(i,:) = F0Electrodes(:,Time)' ;
    else
        %store PVT data for each trial (note cols are pressure, vibration & temperature)
        PVT(i,:) = [ F1pdc(Time)  F1pac(2,Time) F1tdc(Time) ];
        %store electrde data for each trial
        Electrode(i,:) = F1Electrodes(:,Time)' ;
    end
end

%creates a subfolder to save data (otherwise subsequent runs of script will present errors)
mkdir Data 
% save PVT and electrode data in .mat files
save Data/F1_PVT PVT
save Data/F1_Elec Electrode

%% 3D scatter plot

clear all;  %clear up workspace

%load stored data
load Data/F1_PVT

numTrials = 10; %number of trials per object

figure;
for i = 1 : numTrials : height(PVT)
        %plot pressure, vibration & temperature measurements for all trials
        scatter3( PVT(i:i+numTrials-1,1) , PVT(i:i+numTrials-1,2) , PVT(i:i+numTrials-1,3) , 'filled' );
        hold on;
        
end

legend( 'Acrylic' , 'Black Foam' , 'Car Sponge' , 'Flour Sack' , 'Kitchen Sponge' , 'Steel Vase' , 'FontSize' , 16 );
xlabel( 'Pressure (AU)' , 'FontSize' , 16 ) 
ylabel( 'Vibrations (AU)' , 'FontSize' , 16 ) 
zlabel( 'Temperature (AU)' , 'FontSize' , 16 )
