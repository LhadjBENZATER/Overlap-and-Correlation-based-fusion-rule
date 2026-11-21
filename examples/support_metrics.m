function m = support_metrics(truth, decision)
% Bin-level scoring, without frequency tolerance or guard bins.
truth=logical(truth(:)); decision=logical(decision(:));
assert(numel(truth)==numel(decision));
tp=sum(truth & decision); fp=sum(~truth & decision); fn=sum(truth & ~decision);
m=[tp/max(sum(truth),1), fp/max(sum(~truth),1), 2*tp/max(2*tp+fp+fn,1)];
end
