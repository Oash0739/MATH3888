function [tstart,tend,avgstart,avgend,fstart,fend,Hyd1startn,Hyd1endn,mxs,mxe,iFXs,iFXe,dstart,dend]= FFTFunction(time,amps,cuthigh)
t = time;
Hyd1 = amps;

%Before drug is introduced
tstart = t(50:50 + 250);
fst = 1/((tstart(end,1)-tstart(1,1))/(length(tstart)-1));
Hyd1start = Hyd1(50:50 + 250);
avgstart = mean(Hyd1start);
Hyd1startn = Hyd1start - (sum(Hyd1start)/length(Hyd1start));
normal = max(abs(Hyd1startn));
Hyd1startn = Hyd1startn/normal;

%After drug is introduced
tend = t(-250 - 50 + end:end -50);
fen = 1/((tend(end,1)-tend(1,1))/(length(tend)-1));
Hyd1end = Hyd1(-250 - 50 + end:end -50);
avgend = mean(Hyd1end);
Hyd1endn = Hyd1end  - (sum(Hyd1end)/length(Hyd1end));
Hyd1endn = Hyd1endn/normal;

nfft = 100000;

Xs = fft(Hyd1startn,nfft);
Xs = fftshift(Xs);
Xe = fft(Hyd1endn,nfft);
Xe = fftshift(Xe);


fstart = linspace(-nfft/2,nfft/2-1,length(Xs))*fst/nfft;
fend = linspace(-nfft/2,nfft/2-1,length(Xe))*fen/nfft;




mxs = abs(Xs);
mxe = abs(Xe);

%

cstart = [];
for i= 1:length(mxs)-1
    if mxs(i+1) > mxs(i)
        v = [fstart(i + 1); mxs(i + 1)];
        cstart = [cstart v];
    end
end

hfstart = cstart(:,cstart(2,:) > cuthigh*max(mxs));
freqsstart = hfstart(1,:);
freqstart = freqsstart(freqsstart>0);
powerstart = hfstart(2,:);
powerstart = powerstart(freqsstart>0);


cstart = [];
for i= 1:length(freqstart)-1
    if freqstart(i+1) > freqstart(i)
        v = [freqstart(i + 1); powerstart(i + 1)];
        cstart = [cstart v];
    end
end
d = [];
for i = 1:length(cstart(1,:))-1
   if abs(cstart(1,i+1)-cstart(1,i))>0.002
       d = [d cstart(:,i)];
   end
end
dstart = [d cstart(:,end)];


frequencystart = dstart(1,end);




cend = [];
for i= 1:length(mxe)-1
    if mxe(i+1) > mxe(i)
        v = [fend(i + 1); mxe(i + 1)];
        cend = [cend v];
    end
end

hfend = cend(:,cend(2,:) > cuthigh*max(mxe));
freqsend = hfend(1,:);
freqend = freqsend(freqsend>0);
powerend = hfend(2,:);
powerend = powerend(freqsend>0);


cend = [];
for i= 1:length(freqend)-1
    if freqend(i+1) > freqend(i)
        v = [freqend(i + 1); powerend(i + 1)];
        cend = [cend v];
    end
end
d = [];
for i = 1:length(cend(1,:))-1
   if abs(cend(1,i+1)-cend(1,i))>0.002
       d = [d cend(:,i)];
   end
end
dend = [d cend(:,end)];


frequencyend = dend(1,end);




Xs(abs(fstart)>1.1*frequencystart)=0;
Xe(abs(fend)>1.1*frequencyend)=0;


iFXs = ifft(ifftshift(Xs),'symmetric');
iFXs = iFXs(1:length(tstart));
iFXe = ifft(ifftshift(Xe),'symmetric');
%RFXe = real(iFXe);
iFXe = iFXe(1:length(tend));


