% Lightweight integration test of the Figure 11 generator, not a physics simulation.
% A deterministic run_droplet stub forces the "always escapes" branch at every point.
test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(test_dir);
assert(exist(fullfile(repo_root, 'Generate_numerics_data_Fig11.m'), 'file') == 2);
tmp = tempname();
mkdir(tmp);
mkdir(fullfile(tmp, 'Data'));
mkdir(fullfile(tmp, 'Auxillary_functions'));
copyfile(fullfile(repo_root, 'Generate_numerics_data_Fig11.m'), ...
         fullfile(tmp, 'Generate_numerics_data_Fig11.m'));
copyfile(fullfile(test_dir, 'fixtures', 'run_droplet.m'), fullfile(tmp, 'run_droplet.m'));
original_dir = pwd();
cd(tmp);
try
    % Explicit synthetic settings; never use this fixture for published predictions.
    V = 0.2; theta_r = 0; N = 8; tend = 0.001; logout = 0;
    run('Generate_numerics_data_Fig11.m');
    output = fullfile('Data', 'numerics_comparison_data_regenerated.mat');
    assert(exist(output, 'file') == 2, 'Missing regenerated test output');
    saved = load(output);
    assert(isequal(size(saved.xu0_escape), [5, 11]), 'Result array was overwritten');
    assert(all(saved.xu0_escape(:) == -1), 'Expected 55 synthetic escapes');
    assert(isequal(size(saved.hysparmax_out), [1, 5]));
    assert(isequal(size(saved.nu_out), [1, 11]));
    assert(exist(fullfile('Data', 'numerics_comparison_data.mat'), 'file') ~= 2, ...
           'Original reference filename must not be overwritten');
    fprintf('Figure 11 synthetic generator regression passed (55/55 cases).\n');
catch err
    cd(original_dir);
    if exist(tmp, 'dir'), rmdir(tmp, 's'); end
    rethrow(err);
end
cd(original_dir);
rmdir(tmp, 's');
