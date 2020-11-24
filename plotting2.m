%n=1->3 Hydroxycholoroquine
%n=4->6 Caffeine
%n=7->9 Heparin 0.1 umg
%n=10->12 Heparin 0.5 umg
%n=13->15 Heparin 1.0 umg
%n=16->18 Heparin 2.0 umg

frequencies = [];
for n = 1:18
    n = 11;
    cuthigh = 0.99; %percentage height of the dominate frequency
    filename = 'data.xlsx';
    T = xlsread(filename,n);
    % T = cell2mat(t);
    [tstart,tend,avgstart,avgend,fstart,fend,DrugStartNorm,DrugEndNorm,FreqAmpsStartNorm,FreqAmpsEndNorm,IFFTStart,IFFTEnd,DomStart,DomEnd] = FFTFunction(T(:,1),T(:,2),cuthigh);
    v = [DomStart ; DomEnd];
    frequencies = [frequencies v];
end
frequencies = [1:18 ; frequencies];

if n>=1 && n<=3
    str = '50 \mu M Hydroxycholoroquine';
elseif n>3 && n<=6
    str = '40 mM Caffeine';
elseif n>6 && n<=9
    str = 'Heparin 0.1 \mu g mL^{-1}';
elseif n>9 && n<=12
    str ='Heparin 0.5\mu g mL^{-1}';
elseif n>12 && n<=15
    str ='Heparin 1\mu g mL^{-1}';
elseif n>15
    str = 'Heparin 2\mu g mL^{-1}';
end
%%
%Pulling in data from MATCONT
if n>=1 && n<=3
    load('dator\before_hydroxycholoroquine.mat');
elseif n>3 && n<=6
    load('dator\before_caffeine.mat');
elseif n>6 && n<=9
    load('dator\before_0.1_heparin.mat');
elseif n>9 && n<=12
    load('dator\before_0.5_heparin.mat');
elseif n>12 && n<=15
    load('dator\before_1.0_heparin.mat');
elseif n>15
    load('dator\before_2.0_heparin.mat');
end
[tstart1,tend1,avgstart1,avgend1,fstart1,fend1,DrugStartNorm1,DrugEndNorm1,FreqAmpsStartNorm1,FreqAmpsEndNorm1,IFFTStart1,IFFTEnd1,DomStart1,DomEnd1] = FFTFunction(t,y(:,3),cuthigh);
%%
%Plotting normal Graph
figure(1)
plot(T(:,1),T(:,2))
xlabel('time(s)')
ylabel('F-380 Intensity  (arb.)')
title(str)
%%
%Plotting From before and after adding the drug Normalised
figure (2)
subplot(2,3,1)
plot(tstart,DrugStartNorm)
xlabel('time(s)')
ylabel('Intensity (arb.)')
k = ['Before ' str ' added'];
title(join(k))

subplot(2,3,4)
plot(tend,DrugEndNorm)
xlabel('time(s)')
ylabel('Intensity (arb.)')
k = ['After ' str ' added'];
title(join(k))

%Frequency Plots
subplot(2,3,2)
plot(fstart,FreqAmpsStartNorm)
xlabel('Frequencies (Hz)')
ylabel('Frequency weighting')
k = ['Before ' str ' added'];
title(join(k))
% for n = 1:length(DomStart)
%     xline(DomStart(1,n),'--r')
%     xline(-DomStart(1,n),'--r')
% end
yline(cuthigh.*max(FreqAmpsStartNorm),'-r','LineWidth', 1)
xlim([-0.5 0.5])

subplot(2,3,5)
plot(fend,FreqAmpsEndNorm)
xlabel('Frequencies (Hz)')
ylabel('Frequency weighting')
k = ['After ' str ' added'];
title(join(k))
% for n = 1:length(DomEnd)
%     xline(DomEnd(1,n),'--r')
%     xline(-DomEnd(1,n),'--r')
% end
yline(cuthigh.*max(FreqAmpsEndNorm),'-r','LineWidth', 1)
xlim([-0.5 0.5])

%Plotting the IFFT of the cleared data
subplot(2,3,3)
plot(tstart,IFFTStart)
xlabel('time(s)')
ylabel('Intensity (arb.)')
k = ['Before ' str ' added'];
title(join(k))

subplot(2,3,6)
plot(tend,IFFTEnd)
xlabel('time(s)')
ylabel('Intensity (arb.)')
k = ['After ' str ' added' ];
title(join(k))

figure(3)
subplot(2,1,1)
plot(tstart,DrugStartNorm,'-r',tstart,IFFTStart,'-b')
xlabel('time(s)')
ylabel('Intensity (arb.)')
k = ['Before ' str ' added' ];
title(join(k))
k1 = [str ' before clearing'];
k2 = [str ' after clearing'];
legend(k1,k2)

subplot(2,1,2)
plot(tend,DrugEndNorm,'-r',tend,IFFTEnd,'-b')
xlabel('time(s)')
ylabel('Intensity (arb.)')
k = ['After ' str ' added' ];
title(join(k))
k1 = [str ' before clearing'];
k2 = [str ' after clearing'];
legend(k1,k2)

figure(4)
plot(fstart,FreqAmpsStartNorm,fend,FreqAmpsEndNorm)
xlabel('Frequencies (Hz)')
ylabel('Frequency weighting')
k1 = ['before ' str ' is added'];
k2 = ['after ' str ' is clearing'];
legend(k1,k2)
% title(join(k))
xlim([-0.1 0.1])
%%
figure(5)
plot(fstart,FreqAmpsStartNorm,fstart1,FreqAmpsStartNorm1)
xlim([-0.5 0.5])
%%
%Plotting MATCONT Data
meany = mean(y(100:end,3));
topy = max(abs(y(:,3)));
topgarbage = max(abs(y(:,3)-meany));
figure(6)
plot(t,(y(:,3)-meany)/topgarbage*topy,'-r',tend-750,IFFTEnd,'-b')
%%
% Percentage error on model to raw data dominate frequency
if DomEnd1(1,1)>DomEnd(1,1)
    Percent_Error  = abs((DomEnd(1,1) - DomEnd1(1,1))/DomEnd(1,1))
else
    Percent_Error  = abs((-DomEnd(1,1) + DomEnd1(1,1))/DomEnd(1,1))
end

