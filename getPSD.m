function [freq,PSD,DegOfFreedomNum] = getPSD(Input,TimeBinSize)
%UNTITLED Summary of this function goes here
% Input -  raw vector array of the input signal
% TimeBinSize - size of the bin in the time units, i.e. 1 s = 1, 1 ms = 0.001 etc.

SignalLength = length(Input); % bin lenght of the signal
SignalCentered=Input-sum(Input)/SignalLength; % centering the signal by substraction of its mean value

MinWindowNum=3; % Minimal number of windows; can be changed

SpectralRatio =3; % the value determines the number of frequency periods covered by one window, 
%this value must exceed one in order to achieve a satisfactory ratio between variance and bias; can be changed

WindowMinSize= 101; % minimal size of the window

%Note, all the calculation are carried out using the frequency in the units of SignalLength integer fractions

WindowMaxSize=floor((SignalLength-MinWindowNum)/(2*MinWindowNum))*2+1; % maximal window size defined by minimal windows number
MinFreq = round(SpectralRatio/WindowMaxSize*SignalLength); % min frequency defined by signal len, max window size and spectral ratio
MaxFreq = round(0.4*SignalLength); % max frequency, 0.4 may be slightly changed
LogFreqRange = log10(MinFreq):0.05:log10(MaxFreq); %creating frequency range in logarithmic scale with 0.05 step in log scale
FreqRange = ceil(10.^LogFreqRange);% return to the linear scale; 
%Note, it is the legacy code so one can just use ceil(logspace(log10(MinFreq),log10(MaxFreq),NumOfFreqPoints))

WindowSizeArray=max(floor(SpectralRatio*SignalLength./(2*FreqRange))*2+1,WindowMinSize); %size of the windows for each frequency
WindowSizeArray(1)=WindowMaxSize; % size of the window for the minimal windows number (maximal window size)

CycleFreq = 2*pi*FreqRange/SignalLength; % create the cycling frequency for futher calculations
freq=FreqRange/(SignalLength*TimeBinSize); % recalculated freq in Hz units

PSD = CycleFreq*0; %initializing PSD value vector
DegOfFreedomNum=PSD; %initializing the vector that corresponds to the number of degrees of freedom for chi-squared distribution
WinNum=floor(2*SignalLength./WindowSizeArray-1);% number of windows for each frequency (they are overlap so *2)

for ii = 1:length(CycleFreq) % for every frequency
    nn=1:WindowSizeArray(ii); %area covered by one particular window
    WinNumFreq=WinNum(ii); %number of the windows for this particular frequency
    ProfileSize=(WindowSizeArray(ii)-1)/2;% size of the window profile, may be changed for another window profile
    WindowProfile=0.54-0.46*cos(pi*(nn-1)/ProfileSize); % Window profile, that can be changed. E.g. the Hamming window.
    NormalizationCoef=sum(WindowProfile.^2); % normalization coefficint of the Welch's method
    EachWindowPSD=1:WinNumFreq; % initializing the vector of the PSD for each particular window
    for k=1:WinNumFreq % calculating PSD for every window
        EachWindowPSD(k)=abs(sum(exp(nn*1i*CycleFreq(ii)).*WindowProfile.*SignalCentered(floor((k-1)*WindowSizeArray(ii)/2)+nn)))^2/NormalizationCoef; 
    end
    MeanWindowPSD=sum(EachWindowPSD)/WinNumFreq; % mean PSD of all windows
    PSD(ii) = TimeBinSize*MeanWindowPSD/(2*pi); % recalculated PSD; units are Hz^-1
    DegOfFreedomNum(ii)=36*WinNumFreq.^2./(19*WinNumFreq-1); % degrees of freedom number. Here is an example for hamming ang hanning windows.
end
end