function pgcsv(varargin)
% PGCSV - A function for converting geochemical template files in excel
% into a format that is easier to parse into PostgreSQL

%To add: A check on latitude/longitude, some may just specify 'see map/map'
%etc (as ive been doing for some), so would be good to track this so we
%know which ones we can expand on later for this.
%datapath = 'D:\Users\Matt\Documents\University\Doctor of Philosophy\Work for derrick\chemdata\spreadsheets\';
%outpath = 'D:\Users\Matt\Documents\University\Doctor of Philosophy\Work for derrick\heat_production\database\toadd\';
%outpath_dup = 'D:\Users\Matt\Documents\University\Doctor of Philosophy\Work for derrick\heat_production\database\added\';

%outpath = 'D:\Users\Matt\Documents\University\Doctor of Philosophy\Work for derrick\heat_production\database\toadd\';


%Maclaptop
%datapath = '/Users/mgard/Documents/University/PhD/Work for derrick/chemdata/spreadsheets/';
%outpath = '/Users/mgard/Documents/University/PhD/Work for derrick/heat_production/database/toadd/';

%Uni mac
%outpath = '~/Documents/University/PhD/Work for derrick/heat_production/database
% outpath_dup = '/Users/a1628451/Documents/PhD/Work for derrick/heat_production/database/added/';

%Just go back two directories (Reverse / for windows - add option)
datapath = '../../../chemdata/to_add/';
%datapath = '../../../chemdata/spreadsheets/';
%datapath = '../../chemdata/to_fix/';

outpath = '../../../chemdata/pgcsv_outfiles/';
%outpath = '../../heat_production/database/toadd/';

required = {'author','title','year','journal','doi','bibtex','lith_name'...
    'lith_alt_name','prov_name'};
 
if nargin == 1
    d.isdir = 0;
    d.name = varargin{1};
else
    d = dir(datapath);
end

%Check for duplicates - dont run if they exist in the added already
%Optional
%d_dup = dir(outpath_dup);


h = waitbar(0,'Reading files...');
for i = 1:length(d)
    if d(i).isdir || (length(d(i).name) < 5)
        continue;
    elseif ~strcmpi(d(i).name(end-2:end),'xls') ...
        && ~strcmpi(d(i).name(end-3:end),'xlsx')
        waitbar(i/length(d),h,['Skipping ',d(i).name,'...']);
        continue;
    end
    
%Check for duplicates - dont run if they exist in the added already
%Optional
%     flag = 0;
%     if strcmpi(d(i).name(end-2:end),'xls')
%         for k = 3:length(d_dup)
%             if strcmpi(d(i).name(1:end-4),d_dup(k).name(1:end-4))
%                 flag = 1;
%                 break;
%             end
%         end
%     elseif strcmpi(d(i).name(end-3:end),'xlsx')
%         for k = 3:length(d_dup)
%             if strcmpi(d(i).name(1:end-5),d_dup(k).name(1:end-4))
%                 flag = 1;
%                 break;
%             end
%         end
%     end
%     
%     if flag == 1;
%         waitbar(i/length(d),h,['Skipping ',d(i).name,'...']);
%         continue;
%     end
    
    
    
    % load xls(x)
    waitbar((i-0.75)/length(d),h,['Reading ',d(i).name,'...']);
    [num,txt,raw] = xlsread([datapath,d(i).name]);
    [nrt,nct] = size(txt);
    [nrn,ncn] = size(num);

    if nrt ~= nrn
        num = [nan([nrt-nrn,ncn]); num];
    end
    
    %Remove all strange column names
    %e.g. replace \ / with _. Maybe more?
    
    %Get the row names
    temp = txt(22:end,1);
    
    %Check for / or \ and replace.
    codedstring = '\u0338';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,'_');

    codedstring = '\u002F';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,'_');

    
    %Subscript 2
    codedstring = '\u2082';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,'2');
    
    %Subscript 3
    codedstring = '\u2083';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,'3');
    
    %Subscript 5
    codedstring = '\u2085';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,'5');
    
    %Linebreaks etc
    codedstring = '\u2028';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,' ');
    codedstring = '\u000b';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,' ');
    codedstring = '\u000d';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,' ');
    codedstring = '\u0085';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,' ');
    codedstring = '\u2029';
    decodedstring = sprintf(strrep(codedstring, '\u', '\x'));
    temp = regexprep(temp,decodedstring,' ');
    
    %All other symbols
    temp = regexprep(temp,'[^a-zA-Z0-9 _]','');
    
    txt(22:end,1) = temp;
    

    % parse xls(x), return data
    fprintf('%s\n',d(i).name)
    waitbar((i-0.5)/length(d),h,['Parsing ',d(i).name,'...']);
    data = parsexls(num,txt,raw,required);
    
    % write data to csv file
    waitbar((i-0.25)/length(d),h,['Writing ',d(i).name,'...']);
    if strcmpi(d(i).name(end:end-2),'xls')
        fname = [outpath,d(i).name(1:end-3),'csv'];
    else
        fname = [outpath,d(i).name(1:end-4),'csv'];
    end

    % Country replace - replace all country names with 2 letter codes for
    % consistency
    %data.country = country_code(data.country);
    
    writetable(data,fname);

    % fix filename for unix commands
    fname2 = fname;
    ind = strfind(fname2,' ');
    for j = 1:length(ind)
        fname2 = [fname2(1:ind(j)-1),'\',fname2(ind(j):end)];
    end
    
    ind = strfind(fname2,'&');
    for j = 1:length(ind)
        fname2 = [fname2(1:ind(j)-1),'\',fname2(ind(j):end)];
    end
    %NaN's are fine to be included now, I'm dealing with it in python so
    %its universal and doesnt require unix code
    %unix(['sed -i '''' -e ''s/NaN//g'' ',fname2]);
end

close(h);

return


function data = parsexls(num,txt,raw,required)


timefield = {'age','age_min','age_max','age_sd','mage','mage_min','mage_max','mage_sd'};
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
    if ~strcmpi(txt{i,1},'reference') && c == 0
        continue;
    elseif isempty(txt{i,1})
        continue;
    else
        c = c + 1;
    end
    colname{c} = lower(txt{i,1});
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
    
    % fix errors in column names and 
    % rename columns same as postgres database
    if sum(strcmpi(colname{c},{'reference'})) > 0
        colname{c} = 'bibtex';
    end
    
    if sum(strcmpi(colname{c},{'or'})) > 0
        %Cannot have a column named 'or' - conflict with IF/OR etc.
        colname{c} = 'orth';
    end
    
    if sum(strcmpi(colname{c},{'group'})) > 0
        %Cannot have a column named 'group' - conflict with postgresql
        colname{c} = 'chem_group';
    end
    
    if sum(strcmpi(colname{c},{'sfeo'})) > 0
        colname{c} = 'feo_tot';
    end
    
    if sum(strcmpi(colname{c},{'feot'})) > 0
        colname{c} = 'feo_tot';
    end
    
    if sum(strcmpi(colname{c},{'feo_t'})) > 0
        colname{c} = 'feo_tot';
    end
    
    if sum(strcmpi(colname{c},{'fe2o3t'})) > 0
        colname{c} = 'fe2o3_tot';
    end
    
    if sum(strcmpi(colname{c},{'fe2o3_t'})) > 0
        colname{c} = 'fe2o3_tot';
    end
    
    if sum(strcmpi(colname{c},{'feo*'})) > 0
        colname{c} = 'feo_tot';
    end
    
    if sum(strcmpi(colname{c},{'bibkey'})) > 0
        colname{c} = 'bibtex';
    end
    
    if sum(strcmpi(colname{c},{'h2o+'})) > 0
        colname{c} = 'h2o_plus';
    end
    
    if sum(strcmpi(colname{c},{'h2o-'})) > 0
        colname{c} = 'h2o_minus';
    end
    
    if sum(strcmpi(colname{c},{'h2o?'})) > 0
        colname{c} = 'h2o_minus';
    end

    if sum(strcmpi(colname{c},{'geologic_province'})) > 0
        colname{c} = 'prov_name';
    end

    if sum(strcmpi(colname{c},{'subprovince_or_terrane'})) > 0
        colname{c} = 'subprov1';
    end

    if sum(strcmpi(colname{c},{'lithologic_unit'})) > 0
        colname{c} = 'lith_name';
    end
    
    if sum(strcmpi(colname{c},{'rock_type'})) > 0
        colname{c} = 'rock_name';
    end
    
    if sum(strcmpi(colname{c},{'rock type'})) > 0
        colname{c} = 'rock_name';
    end
    
    if sum(strcmpi(colname{c},{'rock facies'})) > 0
        colname{c} = 'rock_composition';
    end
    
    if sum(strcmpi(colname{c},{'rock_facies'})) > 0
        colname{c} = 'rock_composition';
    end

    if sum(strcmpi(colname{c},{'mage_ma'})) > 0
        colname{c} = 'mage';
    end

    if sum(strcmpi(colname{c},{'mage_min_ma'})) > 0
        colname{c} = 'mage_min';
    end

    if sum(strcmpi(colname{c},{'mage_max_ma'})) > 0
        colname{c} = 'mage_max';
    end

    if sum(strcmpi(colname{c},{'mage_sd_ma'})) > 0
        colname{c} = 'mage_sd';
    end

    if sum(strcmpi(colname{c},{'mgo#','mg#','mg_#','Mg_number','mg number','mg #','#mg'})) > 0
        colname{c} = 'mg_num';
    end

    if strcmpi(colname{c},{'sum'})
        colname{c} = 'total';
    end

    if strcmpi(colname{c},{'l.o.i.'})
        colname{c} = 'loi';
    end
    
    if strcmpi(colname{c},{'l.o.i'})
        colname{c} = 'loi';
    end

    if strcmpi(colname{c},{'sample_id'})
        colname{c} = 'sample_name';
    end

    if sum(strcmpi(colname{c},{'na20','2O','20','2o'})) > 0
        colname{c} = 'na2o';
    end

    if strcmpi(colname{c},{'k20'})
        colname{c} = 'k2o';
    end

    if strcmpi(colname{c},{'a1203','al203'})
        colname{c} = 'al2o3';
    end

    if strcmpi(colname{c},{'si02'})
        colname{c} = 'sio2';
    end

    if strcmpi(colname{c},{'ti02'})
        colname{c} = 'tio2';
    end

    if strcmpi(colname{c},{'p205'})
        colname{c} = 'p2o5';
    end

    if strcmpi(colname{c},{'ca0'})
        colname{c} = 'cao';
    end

    if strcmpi(colname{c},{'mg0'})
        colname{c} = 'mgo';
    end

    if strcmpi(colname{c},{'mn0'})
        colname{c} = 'mno';
    end

    if strcmpi(colname{c},{'fe0'})
        colname{c} = 'feo';
    end

    if strcmpi(colname{c},{'fe203'})
        colname{c} = 'fe2o3';
    end
    
    if strcmpi(colname{c},{'Fo203'})
        colname{c} = 'fe2o3';
    end
    
    if strcmpi(colname{c},{'87Rb_86Sr'})
        colname{c} = 'rb87_sr86';
    end
    if strcmpi(colname{c},{'87Rb_86Sr_uncertainty'})
        colname{c} = 'rb87_sr86_uncertainty';
    end
    
    if strcmpi(colname{c},{'87Sr_86Sr'})
        colname{c} = 'sr87_sr86';
    end
    if strcmpi(colname{c},{'87Sr_86Sr_uncertainty'})
        colname{c} = 'sr87_sr86_uncertainty';
    end
    
    if strcmpi(colname{c},{'147Sm_144Nd'})
        colname{c} = 'sm147_nd144';
    end
    if strcmpi(colname{c},{'147Sm_144Nd_uncertainty'})
        colname{c} = 'sm147_nd144_uncertainty';
    end
    
    if strcmpi(colname{c},{'143Nd_144Nd'})
        colname{c} = 'nd143_nd144';
    end
    if strcmpi(colname{c},{'143Nd_144Nd_uncertainty'})
        colname{c} = 'nd143_nd144_uncertainty';
    end
    
    if strcmpi(colname{c},{'176Hf_177Hf'})
        colname{c} = 'hf176_hf177';
    end
    if strcmpi(colname{c},{'176Hf_177Hf_uncertainty'})
        colname{c} = 'hf176_hf177_uncertainty';
    end
    
    if strcmpi(colname{c},{'176Lu_177Hf'})
        colname{c} = 'lu176_hf177';
    end
    if strcmpi(colname{c},{'176Lu_177Hf_uncertainty'})
        colname{c} = 'lu176_hf177_uncertainty';
    end
    
    if strcmpi(colname{c},{'206Pb_204Pb'})
        colname{c} = 'pb206_pb204';
    end
    if strcmpi(colname{c},{'206Pb_204Pb_uncertainty'})
        colname{c} = 'pb206_pb204_uncertainty';
    end
    
    if strcmpi(colname{c},{'207Pb_204Pb'})
        colname{c} = 'pb207_pb204';
    end
    if strcmpi(colname{c},{'207Pb_204Pb_uncertainty'})
        colname{c} = 'pb207_pb204_uncertainty';
    end
    
    if strcmpi(colname{c},{'208Pb_204Pb'})
        colname{c} = 'pb208_pb204';
    end
    if strcmpi(colname{c},{'208Pb_204Pb_uncertainty'})
        colname{c} = 'pb208_pb204_uncertainty';
    end
    
    if strcmpi(colname{c},{'238U_232Th'})
        colname{c} = 'u238_th232';
    end
    if strcmpi(colname{c},{'238U_232Th_uncertainty'})
        colname{c} = 'u238_th232_uncertainty';
    end
    
    if strcmpi(colname{c},{'230Th_232Th'})
        colname{c} = 'th230_th232';
    end
    if strcmpi(colname{c},{'230Th_232Th_uncertainty'})
        colname{c} = 'th230_th232_uncertainty';
    end
    
    if strcmpi(colname{c},{'230Th_238U'})
        colname{c} = 'th230_u238';
    end
    if strcmpi(colname{c},{'230Th_238U_uncertainty'})
        colname{c} = 'th230_u238_uncertainty';
    end
    
    if strcmpi(colname{c},{'226Ra_230Th'})
        colname{c} = 'ra226_th230';
    end
    if strcmpi(colname{c},{'226Ra_230Th_uncertainty'})
        colname{c} = 'ra226_th230_uncertainty';
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
    elseif strfind(colname{c},'_per') | strfind(colname{c},'_wtper')
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
            try
            if ~isnan(num(i,j))
                tmp{j,c} = num(i,j);
            else
                tmp{j,c} = txt{i,j+1};
            end
            catch
                size(num)
                nc,i,j
                error('failed');
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
