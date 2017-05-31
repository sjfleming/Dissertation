%%

load('/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/20151015/20151015_27A_40mV.mat')

%%

delT = 20; % in ms
maxT = 1000; % in ms

tau = 57.81; % fit parameter from the logy, linx fit

xlin = delT/2:delT:maxT;
yy = hist(escapeTimes,xlin);

figure(1)
clf
errorbar(xlin,yy,sqrt(max(1,yy)),'ok')
set(gca,'fontsize',18,'OuterPosition',[0.01 0.01 0.98 0.98],'LooseInset',[0 0 0 0])
xlim([0 maxT])
ylim([0 ceil(max(yy+sqrt(yy))/50)*50])
xlabel('Escape time (ms)')
ylabel('Number of molecules')
title('Histogram of poly(dA) escape times at 60mV, 25°C')
set(gcf,'position',[-878   793   618   297])

%%

%ft = fittype('a * exp(-x/b)');
%f = fit(xlin',yy',ft,'weights',1./max(1,sqrt(yy)),'startpoint',[max(yy), 55]);
ft = fittype('a * exp(-x/57.81)');
f = fit(xlin',yy',ft,'weights',1./max(1,sqrt(yy)),'startpoint',max(yy));
xx = linspace(0,maxT,1000);
figure(1)
hold on
plot(xx,f(xx),'r-')

%%

%yylog = log10(hist(escapeTimes,xlin));
yyn = yy/sum(yy);
lowers = yyn-max(yyn-sqrt(yy)/sum(yy),1e-5);

figure(2)
clf
%errorbar(xlin(~isnan(yylog)),yylog(~isnan(yylog)), ...
%    yylog(~isnan(yylog))-log10(yy(~isnan(yylog))-sqrt(yy(~isnan(yylog)))),log10(yy(~isnan(yylog))+sqrt(yy(~isnan(yylog))))-yylog(~isnan(yylog)),'ok')
errorbar(xlin,yyn,lowers,sqrt(yy)/sum(yy),'ok')
set(gca,'fontsize',18,'OuterPosition',[0.01 0.01 0.98 0.98],'LooseInset',[0 0 0 0],'yscale','log')
xlim([0 maxT])
xlabel('Escape time (ms)')
ylabel('Fraction of molecules')
title('Histogram of poly(dA) escape times at 60mV, 25°C')
set(gcf,'position',[-878   793   618   297])

%%

% ft = fittype('a * exp(-x/b)');
% f2 = fit(xlin',yyn',ft,'weights',1./max(1,sqrt(yy)),'startpoint',[max(yyn), 55]);
ft = fittype('a * exp(-x/57.81)');
f2 = fit(xlin',yyn',ft,'weights',1./max(1,sqrt(yy)),'startpoint',max(yyn));
xx = linspace(0,maxT,1000);
figure(2)
hold on
plot(xx,f2(xx),'r-')
ylim([min(yyn(yyn~=min(yyn)))/1.1 1])

%%

xedges = logspace(log10(5),log10(maxT),round(maxT/delT));
xsizes = diff(xedges);
xlog = (xedges(1:end-1)+xedges(2:end))/2;
yylog = histcounts(escapeTimes,xedges);
yynlog = yylog./(sum(yylog)*xsizes);
lowers = yylog-max(yylog-sqrt(yylog)./(sum(yylog)*xsizes),1e-6);

figure(3)
clf
errorbar(xlog,yynlog,lowers,sqrt(yylog)./(sum(yylog)*xsizes),'ok')
set(gca,'fontsize',18,'OuterPosition',[0.01 0.01 0.98 0.98],'LooseInset',[0 0 0 0],'yscale','log','xscale','log')
xlim([min(xlog)/1.1 maxT])
xlabel('Escape time (ms)')
ylabel('Probability density')
title('Histogram of poly(dA) escape times at 60mV, 25°C')
set(gcf,'position',[-878   793   618   297])

%%

xx = logspace(log10(5),log10(maxT),1000);
ft = fittype('1/a * exp(-x/57.81)');
howfar = 49;
f3 = fit(xlog(1:howfar)',yynlog(1:howfar)',ft,'weights',1./max(1,sqrt(yylog(1:howfar))),'startpoint',55);
figure(3)
hold on
plot(xx,f3(xx),'r-')
ylim([1e-6 1e-1])

%%

figure(3)
xx2 = logspace(log10(5),log10(xlog(end-14)),500);
area(xx2,f3(xx2))
text(30,1e-5,'90% of events','color','r','fontsize',20)
