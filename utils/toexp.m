function name = toexp(seqname)
% returns a folder, which contains .seq for scan, and correspoding traj and .m file

if nargin<1 || isempty(seqname)
    seqname = 'external';
end

d = folder_to_save_seq();

name = fullfile(d,[seqname,'.seq']);
name = strrep(name,'\','/');
end