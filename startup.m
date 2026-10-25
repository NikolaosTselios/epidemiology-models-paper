%% startup.m
% Αυτό το script προσθέτει όλους τους υποφακέλους του project στο MATLAB path
% και ενεργοποιεί το περιβάλλον για εκτέλεση.

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(project_root));

fprintf('Όλοι οι φάκελοι ενεργοποιήθηκαν: %s\n', project_root);