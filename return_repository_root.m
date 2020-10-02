function rootDir = return_repository_root()
%return_repository_root.m
% Returns the path to the local copy of the MATLAB repository.
% Should work for both, pc- or unix-based operating systems.
%
% Input
% none
%
% Output
% unix/pc path to local copy of MATLAB repository
%
%
% Mikus Grasis, CRL, November 2017

thisFile = mfilename('fullpath');
currDir = fileparts(thisFile);

%rootDir = currDir(1:end-6);
rootDir = currDir;

end