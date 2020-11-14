function [filelist]=getALLfile(path,filelist)
if isfile(path)
    filename=split(path,'\');
    filename=filename(size(filename,1)-1,:);
    filename=split(filename,'.');
    filename=filename(1,1);
    filelist=[filelist;filename];
elseif(isdir(path))
    File = dir(fullfile(path));
    Dirname={File.name}';
    for i=1:length(Dirname)
        if(isequal(Dirname{i,:},'.') ||isequal(Dirname{i,:},'..'))
            continue
        else
            path_new=[path,Dirname{i,1},'\'];
            filelist=getALLfile(path_new,filelist);
        end
    end
end
end