function write_skope_params(seqname,msg)
[res_dir,res_name] = fileparts(toexp(seqname));
res_fullname = fullfile(res_dir,[res_name,'.txt']);
str = [seqname,newline,msg];
fid = fopen(res_fullname, 'w');
fprintf(fid, '%s', str);
fclose(fid); 
end