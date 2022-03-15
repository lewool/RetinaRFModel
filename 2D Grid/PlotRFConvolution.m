%% IMPORT DATA
ImportCell=Data.Ecc5mm.Cell1;
CellLAmpData=ImportCell.ResponseFunctions.LIsoResponse.Amplitude;
CellMAmpData=ImportCell.ResponseFunctions.MIsoResponse.Amplitude;
CellLPhaseData=ImportCell.ResponseFunctions.LIsoResponse.Phase;
CellMPhaseData=ImportCell.ResponseFunctions.MIsoResponse.Phase;
CellLMSumAmpData=ImportCell.ResponseFunctions.LMSumResponse.Amplitude;
CellLMSumPhaseData=ImportCell.ResponseFunctions.LMSumResponse.Phase;
CellLMDiffAmpData=ImportCell.ResponseFunctions.LMDiffResponse.Amplitude;
CellLMDiffPhaseData=ImportCell.ResponseFunctions.LMDiffResponse.Phase;
CellCenterWeight=ImportCell.CenterWeight;
CellSurroundWeight=ImportCell.SurroundWeight;
Freqs=ImportCell.Convolution.SpatialFrequencies;

%CD to plot directory
PlotDirectory=strcat('/Users/lauren/Documents/SCHOOL/Dacey/RF Model/2D Grid/Analysis/Cell Plots/',CellName);
if (exist(PlotDirectory) == 0)
    mkdir(PlotDirectory);
    cd(PlotDirectory);
else
    cd(PlotDirectory);
end

%% PLOT DATA

%CONVOLUTION
SFLabels={'SF=1/128' 'SF=1/64' 'SF=1/32' 'SF=1/16' 'SF=1/8' 'SF=1/4' 'SF=1/2' 'SF=1' 'SF=2' 'SF=4' 'SF=8' 'SF=16' 'SF=32' 'SF=64' 'SF=128'};
centerstrength=strcat('Lc=',num2str(CellCenterWeight(1)),', Mc=',num2str(CellCenterWeight(2)),', Sc=',num2str(CellCenterWeight(3)));
surroundstrength=strcat('Ls=',num2str(CellSurroundWeight(1)),', Ms=',num2str(CellSurroundWeight(2)),', Ss=',num2str(CellSurroundWeight(3)));

figure('Name',strcat(CellName,' Response Convolution'));
subplot(2,3,1);
Ymax=max([CellLAmpData CellMAmpData]);
if Ymax==0
    Ymax=10;
end
PlotLAmp=loglog(Freqs,CellLAmpData,'ro');
hold on;
PlotMAmp=loglog(Freqs,CellMAmpData,'go');
%PlotSAmp=loglog(Freqs,TotalSData.Amplitude,'bo');
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Response amplitude (a.u.)'); 
axis([0.005 150 1 100]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
legend([PlotLAmp,PlotMAmp],'Location','southwest')
set(PlotLAmp,...
    'DisplayName','L',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[204 0 0]/255)
set(PlotMAmp,...
    'DisplayName','M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[119 172 48]/255)
%     set(PlotSAmp,...
%         'DisplayName','S',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[0 0.45 0.74])

subplot(2,3,4);
PlotLPhase=semilogx(Freqs,abs(CellLPhaseData),'ro');
hold on;
PlotMPhase=semilogx(Freqs,abs(CellMPhaseData),'go');
%PlotSPhase=semilogx(Freqs,abs(TotalSData.Phase),'go');
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Phase (degrees)'); 
axis([0.005 150 -20 200]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%legend([PlotLPhase,PlotMPhase],'Location','northeast')
set(PlotLPhase,...
    'DisplayName','L',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[204 0 0]/255)
set(PlotMPhase,...
    'DisplayName','M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[119 172 48]/255)
%     set(PlotSPhase,...
%         'DisplayName','S',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[0 0.45 0.74])

subplot(2,3,2);
Ymax=max([CellLMSumAmpData]);
if Ymax==0
    Ymax=10;
end
PlotTotAmp=loglog(Freqs,CellLMSumAmpData,'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Response amplitude (a.u.)'); 
axis([0.005 150 1 100]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
legend(PlotTotAmp,'Location','southwest')
set(PlotTotAmp,...
    'DisplayName','L+M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.8 .8 .8])

subplot(2,3,5);
PlotTotPhase=semilogx(Freqs,abs(CellLMSumPhaseData),'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Phase (degrees)'); 
axis([0.005 150 -20 200]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%legend(PlotTotPhase,'Location','northeast')
set(PlotTotPhase,...
    'DisplayName','L+M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.8 .8 .8])

subplot(2,3,3);
Ymax=max([CellLMDiffAmpData]);
if Ymax==0
    Ymax=10;
end
PlotLMDiffAmp=loglog(Freqs,CellLMDiffAmpData,'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Response amplitude (a.u.)'); 
axis([0.005 150 1 100]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
legend(PlotLMDiffAmp,'Location','southwest')
set(PlotLMDiffAmp,...
    'DisplayName','L-M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.31 .31 .31])

subplot(2,3,6);
PlotLMDiffPhase=semilogx(Freqs,abs(CellLMDiffPhaseData),'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Phase (degrees)'); 
axis([0.005 150 -20 200]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%legend(PlotTotPhase,'Location','northeast')

set(PlotLMDiffPhase,...
    'DisplayName','L-M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.31 .31 .31])

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,{centerstrength;surroundstrength},'HorizontalAlignment','center','VerticalAlignment', 'top');

saveas(gcf,strcat(CellName,' Response Convolution'));
