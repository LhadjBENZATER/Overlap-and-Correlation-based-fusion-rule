function run_validation
root=fileparts(mfilename('fullpath'));cd(root);addpath(root,fullfile(root,'tests'));
if ~exist('docs','dir'),mkdir('docs');end
diary(fullfile('docs','validation.txt'));
cleanup=onCleanup(@()diary('off'));
fprintf('Validation date (UTC runner): %s\nMATLAB: %s\n',datestr(now,31),version);
fprintf('Scope: numerical checks and synthetic onboarding demo; legacy experiments excluded.\n');
run_tests;
timer=tic;result=demo_fusion(true);
assert(all(isfinite(result.metrics(:))));
fprintf('Demo elapsed seconds: %.6g\n',toc(timer));
fprintf('Validation completed successfully.\n');
end
