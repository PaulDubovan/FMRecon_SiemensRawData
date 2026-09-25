function copymfiletoexp(newname)
% copy the current m-file to a folder for future reference (e.g. for debug,
% to see what changed between different sessions)


fname = matlab.desktop.editor.getActiveFilename;
d = folder_to_save_seq();
newfullname = fullfile(d,[newname,'.m']);
status = copyfile(fname,newfullname);
if status==0
    warning('Failed to copy file from \n%s to\n %s',fname,newfullname);
else
    fprintf('\nFor future reference\n%s \nis copied to \n%s\n', fname,newfullname);
end
end