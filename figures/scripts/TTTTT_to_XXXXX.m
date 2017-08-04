% script to plot TTTTT to XXXXX level changes at different voltages

%% data locations

files = {'C:\Stephen\Research\Data\Biopore\20170727\2017_07_27_0008.abf', ...
    'C:\Stephen\Research\Data\Biopore\20170727\2017_07_27_0005.abf', ...
    'C:\Stephen\Research\Data\Biopore\20170727\2017_07_27_0006.abf', ...
    'C:\Stephen\Research\Data\Biopore\20170727\2017_07_27_0007.abf'};

tr = {[134.3557      138.0492], ...
    [771.3389      773.0883], ...
    [620.3198      623.8075], ...
    [433.251      437.2014]};

V = [80, 100, 140, 160];
I0 = [159.19, 196.96, 275.23, 318.60];

%% grab data

d = cell(numel(V));
for i = 1:numel(V)
    
    sigdata = SignalData(files{i});
    filtname = sprintf('Low-pass Bessel (%d Hz)', 1000);
    fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,1000),filtname);
    fsigs = fsigs(1);
    data = sigdata.getByTime(tr{i});
    d{i} = data(:,[1,fsigs]);
    
    figure(i)
    plot((d{i}(:,1) - d{i}(1,1)),d{i}(:,2)*1000/I0(i))
    xlabel('Time (s)')
    ylabel('Current (I/I_0)')
    title([num2str(V(i)) ' mV'])
    set(gca,'fontsize',14,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
    ylim([0.05 0.45])
    set(gcf,'position',[-948   533   461   420])
    
    figure(i+numel(V))
    digitization = 1/10 * 0.30517578125 * 10; % in pA
    current = d{i}(:,2)*1000; % in pA
    xx = 0:digitization:max(current);
    yy = hist(current,xx);
    barh(xx/I0(i),yy)
    xlabel('Counts')
    ylabel('Current (I/I_0)')
    title([num2str(V(i)) ' mV'])
    set(gca,'fontsize',14,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
    ylim([0.05 0.45])
    set(gcf,'position',[-446   533   340   420])
    
end

%% plot

