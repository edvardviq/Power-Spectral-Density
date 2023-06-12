function [] = plotPSD(freq,PSD,NumOfDeg,color)
%Function to plot the PSD
%upper and lower bound of 95%-confidence interval for chi-squared distribution
Sup=NumOfDeg.*PSD./chi2inv(0.025,NumOfDeg);
Sdn=NumOfDeg.*PSD./chi2inv(0.975,NumOfDeg);
    
errorbar(freq,PSD,PSD-Sdn,Sup-PSD,'s','Color',color,'MarkerSize',10,'LineWidth',1.2);
set(gca,'YScale','log','XScale','log','XLim',[min(freq)/2 max(freq)*2],'YLim',[0.5*min(Sdn),2*max(Sup)])
xlabel('Frequency [Hz]')
ylabel('PSD [Hz^{-1}]')
grid on