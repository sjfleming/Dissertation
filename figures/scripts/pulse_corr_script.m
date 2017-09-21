%%

% open a file in PoreView
% filter data to 500Hz
% select event with cursors
% run this script
% get kappa, the correlation metric

tr = pv.getCursors();
filter = 500;
current = util.downsample_pointwise(pv.data, pv.psigs(1).sigs(end), tr, min(filter*5*diff(tr),10e6));
t = linspace(tr(1),tr(2),numel(current));
data = [t', current'];
clear t current
levels = karplus_levels(data, 1/(1600/1000), 1e-10, filter);
clear data
[k, pulses] = pulse_correlation(pv.data,tr,4,levels);
disp(k)


%%

% 5ms in front and in back, compare

backup = 2e-4;
skip = 2e-3;
comptime = 2e-3;

pulses = [];
tr = pv.getCursors();
d = pv.data.getViewData(tr);
logic = diff(medfilt1(d(:,4),2))>0.5;
logic = logic & ~[false; logic(1:end-1)]; % one point every spike
inds = find(logic); % indices of those spikes
inds = inds([true; diff(inds)>10]);
dt = d(2,1)-d(1,1);
for i = 1:numel(inds)
    ind = (tr(1) + inds(i)*dt)/pv.data.si; % index of raw data
    pulses(end+1) = pv.data.si * pv.data.findNext(@(d) d(:,4)>0.05,ind-100);
end

h = nan(1,numel(pulses));
p = nan(1,numel(pulses));
diff_mean = nan(1,numel(pulses));
val = nan(1,numel(pulses));
for i = 1:numel(pulses)
    d1 = pv.data.getByTime([pulses(i)-backup-comptime, pulses(i)-backup]);
    d2 = pv.data.getByTime([pulses(i)+skip, pulses(i)+skip+comptime]);
    [h(i),p(i)] = kstest2(d1(:,2),d2(:,2));
    diff_mean(i) = mean(d1(:,2))-mean(d2(:,2));
    val(i) = std(d1(:,2)) * std(d2(:,2)) / std([d1(:,2); d2(:,2)]);
end

disp(sum(log10(val)<-3)/numel(pulses))

