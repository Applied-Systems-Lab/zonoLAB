for j = 0:9
    filePath = 'C:\Users\huang\Downloads\mnist_png-master\mnist_png\testing\' + string(j);
    filelist=dir(fullfile(filePath, '*.png'));
    for i = 1:numel(filelist)
        oldname = filelist(i).name;
        newname = [num2str(i) '_1.png']; % can also do newname = string(i) + ".png";
        
        fprintf('will rename %s to %s\n'); %check these before running with the next line uncommented
        movefile(fullfile(filePath, oldname), fullfile(filePath, newname));
    end
end