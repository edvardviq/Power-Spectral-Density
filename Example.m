%Here is the example
% blinking_data=load(char('Example.csv'),'-ascii'); % one may use .csv to load the data
blinking_data=load('Example.txt'); %using .txt exaple blinking trajectory, this trajectory corresponds to the dot 3.2 from the doi 10.1039/D3TC00638G
timeArray=blinking_data(:,1); %time with 0.01 s bin
signalArray=blinking_data(:,2);
[freq,PSD,degrees]=getPSD(signalArray',0.01); %getting the PSD using modified Welch's method described in doi 10.1021/nl3035674 and 10.1002/aenm.202102449
color='k'; %one of standart Matlab colors; one can change it
plotPSD(freq,PSD,degrees,color) %plot the PSD 