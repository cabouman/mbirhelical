% FILENAME:  myscript.m

% Display name of compute node which ran this job.
[c name] = system('hostname');
fprintf('\n\nhostname:%s\n', name);

% Display three random numbers.
A = rand(1,3);
fprintf('%f %f %f\n', A);

quit;
