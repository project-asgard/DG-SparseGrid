function [root_directory] = get_root_folder()

path_to_asgard = which('asgard');

root_directory = fileparts(path_to_asgard);
