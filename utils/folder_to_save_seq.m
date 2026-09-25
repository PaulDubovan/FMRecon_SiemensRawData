function d = folder_to_save_seq()
if ispc || isunix && strcmp(  getenv('USER'),'xy992'  ) % xingwang
    d0 = fullfile('C:\Users\76494\Syncbb', 'seq-for-exp');
    yyyymmdd = datestr(datetime('now'), 'yyyymmdd');
    yyyy = yyyymmdd(1:4);
    d = fullfile(d0, yyyy, yyyymmdd);
elseif isunix && strcmp(  getenv('USER'),'berkin'  ) % berkin
    d = '/cluster/berkin/berkin/Matlab_Code_New/pulseq/spiral/todo';
else % others
    if ispc
        homed = getenv('USERPROFILE');
    else
        homed = getenv('HOME');
    end
    d = fullfile(homed,'todo');
end

if ~exist(d,'dir')
    mkdir(d)
end
end