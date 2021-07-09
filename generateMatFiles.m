errors = [];

for i = 48:size(data, 1)
    try
        [a,~,c] = io.dataFactoryGratingSubspace(i);
    catch
        disp('error in number:')
        disp(i)
        errors = [errors; i];
    end
    close all
    fclose all
end