function slab = loadslab

fid = fopen('../GIS/subduction/ggge761-sup-0006-ts03.txt','r');

str = '>>>>Digitized slab contours for ';
m = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if strcmp(tline(1:3),'>>>')
        m = m + 1;
        slab(m).name = tline(length(str):end);
        c = 1;
        fprintf('Reading slab %s...\n',slab(m).name);
        continue
    elseif strcmp(tline(1:3),'>  ')
        [~,tline] = strtok(tline,' ');
        
        % number of lines to read in contour
        [n,tline] = strtok(tline,' ');
        
        % depth of contour
        [depth,tline] = strtok(tline,' ');
        slab(m).depth(c) = str2num(depth);
    end
    tmp = textscan(fid,'%f%f%f',str2double(n));
    %[tmp{1} tmp{2}]
    slab(m).contour(c).lat = tmp{2};
    slab(m).contour(c).lon = tmp{1};
    c = c + 1;
    
    tline = fgetl(fid);
    if ~ischar(tline), break, end
end
fclose(fid);

return