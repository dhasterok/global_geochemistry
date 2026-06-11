function outfile = proto_xls2csv(varargin)
% PROTO_XLS2CSV - convert data file to format for protolith determination.
%
%   outfile = proto_xls2csv will open a dialog to select an infile '*.xls'
%   or '*.xlsx to convert to a '*.csv' format necessary for protolith class
%   prediction.  The function will output the name and path to the outfile.
%
%   proto_xls2csv(infile) can be used to specify the infile directly instead of
%   using the dialog to open the infile.
%
%   proto_xls2csv(path,infile) can be used to open the infile if it is in a
%   different path than the working directory.
%
%   Note the code will add some additional columns to the outfile that are
%   used by the global geochemistry database processing codes.
%
%   Example:
%       proto_xls2csv('xlsx/','protolith_template.xlsx');

% Original: 2019 June 18 by D. Hasterok dhasterok@gmail.com

% get infile name
switch nargin
    case 0
        [infile,path] = uigetfile({'*.xlsx';'*.xls'},'Select a geochemistry data infile')
    case 1
        infile = varargin{1};
        path = '';
    case 2
        infile = varargin{2};
        path = varargin{1};
        if ~strcmp(path(end),'/')
            path = [path,'/'];
        end
    otherwise
        error('Incorrect number of inputs.');
end
if infile == 0
    error('File not found or invalid infile name.');
end

% outfile name
if strcmpi(infile(end:end-2),'xls')
    outfile = [path,infile(1:end-3),'csv'];
else
    outfile = [path,infile(1:end-4),'csv'];
end

% load infile
fprintf('Loading %s...\n',infile)
[num,txt,raw] = xlsread([path,infile]);
[nrt,nct] = size(txt);
[nrn,ncn] = size(num);

if nrt ~= nrn
    num = [nan([nrt-nrn,ncn]); num];
end

%Get the row names
temp = txt(22:end,1);

%Check for special characters and replace
fprintf('Cleaning %s...\n',infile)
codedstring = {'\u0338', '_'; % /
    '\u002F', '_'; % /
    '\u2080', '0'; % subscript 0
    '\u2081', '1'; % subscript 1
    '\u2082', '2'; % subscript 2
    '\u2083', '3'; % subscript 3
    '\u2084', '4'; % subscript 4
    '\u2085', '5'; % subscript 5
    '\u2086', '6'; % subscript 6
    '\u2087', '7'; % subscript 7
    '\u2088', '8'; % subscript 8
    '\u2089', '9'; % subscript 9
    '\u2028', ' '; % line separator
    '\u000b', ' '; % tab
    '\u000d', ' '; % carriage return
    '\u0085', ' '; % next line
    '\u2029', ' '};% paragraph separator

for i = 1:size(codedstring,1)
    decodedstring = sprintf(strrep(codedstring{i,1}, '\u', '\x'));
    temp = regexprep(temp,decodedstring,codedstring{i,2});
end

%All other symbols
temp = regexprep(temp,'[^a-zA-Z0-9 _]','');
txt(22:end,1) = temp;

% parse xls(x), return data
fprintf('Parsing %s...\n',infile)
data = parsexls(num,txt,raw);

% write data to csv outfile
fprintf('Writing %s...\n',outfile)
writetable(data,outfile);

return


function data = parsexls(num,txt,raw)

% some of this is for the global geochemical database
required = {'author','title','year','journal','doi','bibtex','lith_name'...
    'lith_alt_name','prov_name'};
timefield = {'age','age_min','age_max','age_sd'};
numfield = {'sio2','tio2','al2o3','fe2o3','fe2o3_tot','feo','feo_tot', ...
    'mgo','cao','na2o','k2o','p2o5','mno','sro','bao','cr2o3','nio','caco3','so3', ...
    'h2o_tot','h2o_plus','h2o_minus','co2','loi','total','latitude','longitude'};

% numfield = {'sio2','tio2','al2o3','fe2o3','fe2o3_tot','feo','feo_tot', ...
%     'mgo','cao','na2o','k2o','p2o5','mno','cr2o3','nio','caco3','so3', ...
%     'h2o_tot','co2','loi','total','latitude','longitude'};
% Added: sro, h2o_minus, h2o_plus, 'bao'

el = {'li','na','k','rb','cs', ...
    'be','mg','ca','sr','ba', ...
    'sc','y', ...
    'la','ce','th','pr','nd','u', ...
    'sm','eu','gd','tb','dy','ho', ...
    'er','tm','yb','lu', ...
    'ti','zr','hf','v','nb','ta', ...
    'cr','mo','w','mn','fe','co','ni', ...
    're','ru','os','rh','pd','pt','ir' ...
    'cu','ag','au', ...
    'zn','cd','hg', ...
    'b','al','ga','in','tl', ...
    'c','si','ge','sn','pb', ...
    'n','p','as','sb','bi', ...
    's','sc','se','te', ...
    'f','cl','br','i'};

% Added: sc

[nr,nc] = size(txt);

c = 0;
for i = 1:nr
    if isempty(txt{i,1})
        continue;
    else
        c = c + 1;
    end
    colname{c} = lower(txt{i,1});
    ind = findstr(colname{c},'(%)');
    if ~isempty(ind)
        colname{c} = [colname{c}(1:ind-1),'wtper'];
    end
    ind = findstr(colname{c},'(');
    if ~isempty(ind)
        colname{c} = colname{c}(1:ind-2);
    end
    ind = findstr(colname{c},' ');
    colname{c}(ind) = '_';

    %Remove UTF characters that appear as 'blank' but not removed by
    %strtrim/etc. Move this to somewhere that can deal with descriptions
    %etc as well. I believe in the bdl section but not sure what the 'flag'
    %etc means.
    colname = regexprep(colname,char(13), '');
    colname = regexprep(colname,char(160), '');
    colname = regexprep(colname,char(8199), '');
    colname = regexprep(colname,char(8232), '');
    colname = regexprep(colname,char(8233), '');
    colname = regexprep(colname,'-','');
    colname = regexprep(colname,'—','');

    while 1
        if strcmpi(colname{c}(end),'_')
            colname{c} = colname{c}(1:end-1);
        else
            break;
        end
    end
    
    % fix errors in column names or rename columns same as global
    % geochemical database
    
    % add new column name errors or columns to rename as necessary
    colstring = {'ident', 'sample_name';
        'reference', 'bibtex';
        'or', 'orth';
        'group', 'chem_group';
        'sfeo', 'feo_tot';
        'feot', 'feo_tot';
        'feo_t', 'feo_tot';
        'fe2o3t', 'fe2o3_tot';
        'fe2o3_t', 'fe2o3_tot';
        'feo*', 'feo_tot';
        'bibkey', 'bibtex';
        'h2o+', 'h2o_plus';
        'h2o-', 'h2o_minus';
        'h2o?', 'h2o_minus';
        'geologic_province', 'prov_name';
        'subprovince_or_terrane', 'subprov1';
        'lithologic_unit', 'lith_name';
        'rock_type', 'rock_name';
        'rock type', 'rock_name';
        'rock facies', 'rock_composition';
        'rock_facies', 'rock_composition';
        'mgo#', 'mg_num';
        'mg#', 'mg_num';
        'mg_#', 'mg_num';
        'Mg_number', 'mg_num';
        'mg number', 'mg_num';
        'mg #', 'mg_num';
        '#mg', 'mg_num';
        'sum', 'total';
        'l.o.i.', 'loi';
        'l.o.i', 'loi';
        'sample_id', 'sample_name';
        'na20', 'na2o';
        '2O', 'na2o';
        '20', 'na2o';
        '2o','na2o';
        'k20', 'k2o';
        'a1203', 'al2o3';
        'al203', 'al2o3';
        'si02', 'sio2';
        'ti02', 'tio2';
        'p205', 'p2o5';
        'ca0', 'cao';
        'mg0', 'mgo';
        'mn0', 'mno';
        'fe0', 'feo';
        'fe203', 'fe2o3';
        'Fo203', 'fe2o3';
        '87Rb_86Sr', 'rb87_sr86';
        '87Rb_86Sr_uncertainty', 'rb87_sr86_uncertainty';
        '87Sr_86Sr', 'sr87_sr86';
        '87Sr_86Sr_uncertainty', 'sr87_sr86_uncertainty';
        '147Sm_144Nd', 'sm147_nd144';
        '147Sm_144Nd_uncertainty', 'sm147_nd144_uncertainty';
        '143Nd_144Nd', 'nd143_nd144';
        '143Nd_144Nd_uncertainty', 'nd143_nd144_uncertainty';
        '176Hf_177Hf', 'hf176_hf177';
        '176Hf_177Hf_uncertainty', 'hf176_hf177_uncertainty';
        '176Lu_177Hf', 'lu176_hf177';
        '176Lu_177Hf_uncertainty', 'lu176_hf177_uncertainty';
        '206Pb_204Pb', 'pb206_pb204';
        '206Pb_204Pb_uncertainty', 'pb206_pb204_uncertainty';
        '207Pb_204Pb', 'pb207_pb204';
        '207Pb_204Pb_uncertainty', 'pb207_pb204_uncertainty';
        '208Pb_204Pb', 'pb208_pb204';
        '208Pb_204Pb_uncertainty', 'pb208_pb204_uncertainty';
        '238U_232Th', 'u238_th232';
        '238U_232Th_uncertainty', 'u238_th232_uncertainty';
        '230Th_232Th', 'th230_th232';
        '230Th_232Th_uncertainty', 'th230_th232_uncertainty';
        '230Th_238U', 'th230_u238';
        '230Th_238U_uncertainty', 'th230_u238_uncertainty';
        '226Ra_230Th', 'ra226_th230';
        '226Ra_230Th_uncertainty', 'ra226_th230_uncertainty'};  
        
    for j = 1:size(colstring,1)        
        if strcmpi(colname{c},colstring{j,1})
            colname{c} = colstring{j,2};
        end
    end
    
    
    
    flag = 0;
    %This doesnt capture < symbols in major element columns
    if strfind(colname{c},'_ppm')
        colname{c} = colname{c}(1:end-4);
    elseif strfind(colname{c},'_ppb')
        num(i,:) = num(i,:)*1e3;
        colname{c} = colname{c}(1:end-4);
    elseif strfind(colname{c},'_ppt')
        num(i,:) = num(i,:)*1e-6;
        colname{c} = colname{c}(1:end-4);
    elseif strfind(colname{c},'_per') | strfind(colname{c},'_wtper') ...
            | strfind(colname{c},'_%')
        num(i,:) = num(i,:)*1e4;
        colname{c} = colname{c}(1:end-4);
    end

    if sum(strcmpi(colname{c},el)) > 0
        colname{c} = [colname{c},'_ppm'];
        flag = 1;
    elseif sum(strcmpi(colname{c},numfield)) > 0 ...
        || sum(strcmpi(colname{c},timefield)) > 0
        flag = 1;
    end

    if flag
        for j = 1:length(num(i,:))
            if isnan(num(i,j)) && ~isempty(txt{i,j})
                ind = strfind(txt{i,j},'< ');
                if ~isempty(ind)
                    txt{i,j} = ['-',txt{i,j}(ind+1:end)];

                    [N,ok] = str2num(txt{i,j});
                    if ok
                        num(i,j) = N;
                        continue;
                    end
                end

                % subtle difference to above '< ', which includes a space
                %Consider using strtrim and regexprep?
                ind = findstr(txt{i,j},'<');
                if ~isempty(ind)
                    txt{i,j} = ['-',txt{i,j}(ind+1:end)];

                    [N,ok] = str2num(txt{i,j});
                    if ok
                        num(i,j) = N;
                        continue;
                    end
                end

                % there seem to be so few '>' that we can safely ignore them
                ind = findstr(txt{i,j},'>');
                if ~isempty(ind)
                    txt{i,j}(ind) = '';
                    [N,ok] = str2num(txt{i,j});
                    if ok
                        num(i,j) = N;
                        continue;
                    end
                end

                % if there is a single '-' remove it, note this does not work
                % with the endash.
                if strcmpi(txt{i,j},'-')
                    continue;
                end

                if strcmpi(txt{i,j},'bdl') || strcmpi(txt{i,j},'b.d.l.') ...
                    || strcmpi(txt{i,j},'nd') || strcmpi(txt{i,j},'n.d.')...
                    || strcmpi(txt{i,j},'bd') || strcmpi(txt{i,j},'b.d.')
                    num(i,j) = 0;
                end
            end
        end
        % change < to - , exact '-' to NaN, '>' to blank, 'b.d.l.' and 'n.d.' to 0

        %if sum(strcmpi(colname{c},timefield)) > 0
        %    % create text-based timescale names
        %    colname{c} = 
        %end

        for j = 1:nc-1
            tmp{j,c} = num(i,j);
        end
    else
        for j = 1:nc-1
        %    txtcol{j,1} = txt{i,j+1};
            if ~isnan(num(i,j))
                tmp{j,c} = num(i,j);
            else
                tmp{j,c} = txt{i,j+1};
            end
        end
        %tmp{c} = txtcol;
    end
end

c = 1;
while c < length(colname)
    clear coltmp tmp2
    if sum(strcmpi(colname,colname{c})) == 1
        c = c + 1;
        continue;
    end
    
    ind = find(strcmpi(colname,colname{c}));
    for j = 1:length(ind)
        for i = 1:nc-1
            try
                dup(i,j) = tmp{i,ind(j)};
            catch
                fprintf('Check to see if there is a duplicate text row.\n');
                colname{c}
                [i,ind(j)]
                size(tmp{i,ind(j)})
                tmp{i,ind(j)}
                return
            end
        end
    end

    val = mean(dup,2,'omitnan');

    d = 1;
    for j = 1:length(colname)
        if j == c
            for i = 1:nc-1
                tmp2{i,d} = val(i);
            end
        elseif intersect(j,ind)
            continue;
        else
            for i = 1:nc-1
                tmp2{i,d} = tmp{i,j};
            end
        end

        coltmp{d} = colname{j};
        d = d + 1;
    end
    colname = coltmp;
    tmp = tmp2;

    c = c + 1;
end

[nr,nc] = size(tmp);

x = 1;
for i = 1:length(required)
    if sum(strcmpi(colname,required{i})) == 0
        colname{nc+x} = required{i};
        for i = 1:nr
            tmp{i,nc+x} = '';
        end
        x = x + 1;
    end
end

data = cell2table(tmp,'VariableNames',colname);

return
