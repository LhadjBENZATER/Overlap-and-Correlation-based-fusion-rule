function run_tests
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'examples'));
m=support_metrics([1 1 0 0],[1 0 1 0]);assert(max(abs(m-[0.5 0.5 0.5]))<1e-12);
m=support_metrics([1 1 0 0],[1 1 0 0]);assert(isequal(m,[1 0 1]));
m=support_metrics([1 1 0 0],[0 0 0 0]);assert(isequal(m,[0 0 0]));
assert(isequal(comb([0;1;1;1;0],[0;0;1;0;0],0.5),[0;1;1;1;0]));
assert(isequal(comb([1;1;0],[1;0;0],0.5),[1;1;0]));
assert(isequal(comb([0;1;1],[0;0;1],0.5),[0;1;1]));
assert(~any(comb(ones(4,1),zeros(4,1),0.5)));
fprintf('All numerical checks passed.\n');
end
