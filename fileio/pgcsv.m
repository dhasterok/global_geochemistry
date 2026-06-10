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
    txt = strtrim(txt);

    [nrt,nct] = size(txt)
    [nrn,ncn] = size(num)

    if nrt ~= nrn
        num = [nan([nrt-nrn,ncn]); num];
    end
    
    %Remove all strange column names
    %e.g. replace \ / with _. Maybe more?
    
    %Get the row names
    temp = txt(7:end,1)
    
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
    
    txt(7:end,1) = temp;
    

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
    'h2o_tot','h2o_plus','h2o_minus','co2','loi','total','latitude','longitude','mswd','mmswd'};

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
    % skip rows until the 'reference' is reached as it is the start of the
    % data table.
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
    if any(strcmpi(colname{c},{'reference'}))
        colname{c} = 'bibtex';
    end
    
    if any(strcmpi(colname{c},{'or'}))
        %Cannot have a column named 'or' - conflict with IF/OR etc.
        colname{c} = 'orth';
    end
    
    if any(strcmpi(colname{c},{'group'}))
        %Cannot have a column named 'group' - conflict with postgresql
        colname{c} = 'chem_group';
    end
    
    if any(strcmpi(colname{c},{'sfeo', 'feot', 'feo_t', 'tfeo', 'feo*','feotot'}))
        colname{c} = 'feo_tot';
    end

    if any(strcmpi(colname{c},{'sfe2o3', 'fe2o3t', 'tfe2o3', 'fe2o3_t','fe2o3tot'}))
        colname{c} = 'fe2o3_tot';
    end
    
    if any(strcmpi(colname{c},{'s_tot'}))
        colname{c} = 's_wtper';
    end

    if any(strcmpi(colname{c},{'c_tot'}))
        colname{c} = 'c_wtper';
    end
    
    if any(strcmpi(colname{c},{'bibkey'}))
        colname{c} = 'bibtex';
    end
    
    if any(strcmpi(colname{c},{'h2o+'}))
        colname{c} = 'h2o_plus';
    end
    
    if any(strcmpi(colname{c},{'h2o-'}))
        colname{c} = 'h2o_minus';
    end
    
    if any(strcmpi(colname{c},{'h2o?'}))
        colname{c} = 'h2o_minus';
    end

    if any(strcmpi(colname{c},{'geologic_province'}))
        colname{c} = 'prov_name';
    end

    if any(strcmpi(colname{c},{'subprovince_or_terrane'}))
        colname{c} = 'subprov1';
    end

    if any(strcmpi(colname{c},{'lithologic_unit'}))
        colname{c} = 'lith_name';
    end
    
    if any(strcmpi(colname{c},{'rock_type'}))
        colname{c} = 'rock_name';
    end
    
    if any(strcmpi(colname{c},{'rock type'}))
        colname{c} = 'rock_name';
    end
    
    if any(strcmpi(colname{c},{'rock facies'}))
        colname{c} = 'rock_composition';
    end
    
    if any(strcmpi(colname{c},{'rock_facies'}))
        colname{c} = 'rock_composition';
    end

    if any(strcmpi(colname{c},{'mage_ma'}))
        colname{c} = 'mage';
    end

    if any(strcmpi(colname{c},{'mage_min_ma'}))
        colname{c} = 'mage_min';
    end

    if any(strcmpi(colname{c},{'mage_max_ma'}))
        colname{c} = 'mage_max';
    end

    if any(strcmpi(colname{c},{'mage_sd_ma'}))
        colname{c} = 'mage_sd';
    end

    if any(strcmpi(colname{c},{'mgo#','mg#','mg_#','Mg_number','mg number','mg #','#mg','mg-no'}))
        colname{c} = 'mg_num';
    end

    if strcmpi(colname{c},{'sum'})
        colname{c} = 'total';
    end

    if any(strcmpi(colname{c},{'l.o.i.', 'l.o.i' 'loss'}))
        colname{c} = 'loi';
    end

    if strcmpi(colname{c},{'sample_id'})
        colname{c} = 'sample_name';
    end

    if any(strcmpi(colname{c},{'na20','2O','20','2o'}))
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

    if strcmpi(colname{c},{'84Sr_86Sr'})
        colname{c} = 'sr84_sr86';
    end

    if strcmpi(colname{c},{'84Sr_86Sr_uncertainty'})
        colname{c} = 'sr84_sr86_uncertainty';
    end

    if strcmpi(colname{c},{'84Sr_88Sr'})
        colname{c} = 'sr84_sr88';
    end
    
    if strcmpi(colname{c},{'147Sm_144Nd'})
        colname{c} = 'sm147_nd144';
    end

    if strcmpi(colname{c},{'147Sm_143Nd'})
        colname{c} = 'sm147_nd143';
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

    if strcmpi(colname{c},{'146Nd_145Nd'})
        colname{c} = 'nd146_nd145';
    end
    if strcmpi(colname{c},{'146Nd_145Nd_uncertainty'})
        colname{c} = 'nd146_nd145_uncertainty';
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

    if strcmpi(colname{c},{'187Re_188Os'})
        colname{c} = 're187_os188';
    end
    if strcmpi(colname{c},{'187Re_188Os_uncertainty'})
        colname{c} = 're187_os188_uncertainty';
    end
    
    if strcmpi(colname{c},{'187Os_188Os'})
        colname{c} = 'os187_os188';
    end
    if strcmpi(colname{c},{'187os_188Os_uncertainty'})
        colname{c} = 'os187_os188_uncertainty';
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

    if strcmpi(colname{c},{'206Pb_238U'})
        colname{c} = 'pb206_u238';
    end
    if strcmpi(colname{c},{'206Pb_238U_uncertainty'})
        colname{c} = 'pb206_u238_uncertainty';
    end
    
    if strcmpi(colname{c},{'207Pb_238U'})
        colname{c} = 'pb207_u238';
    end
    if strcmpi(colname{c},{'207Pb_238U_uncertainty'})
        colname{c} = 'pb207_u238_uncertainty';
    end

    if strcmpi(colname{c},{'206Pb_235U'})
        colname{c} = 'pb206_u235';
    end
    if strcmpi(colname{c},{'206Pb_235U_uncertainty'})
        colname{c} = 'pb206_u235_uncertainty';
    end
    
    if strcmpi(colname{c},{'207Pb_235U'})
        colname{c} = 'pb207_u235';
    end
    if strcmpi(colname{c},{'207Pb_235U_uncertainty'})
        colname{c} = 'pb207_u235_uncertainty';
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
    if strcmpi(colname{c},{'238U_204Pb'})
        colname{c} = 'u238_pb204';
    end
    if strcmpi(colname{c},{'232Th_204Pb'})
        colname{c} = 'th232_pb204';
    end
    flag = 0;
    %This doesnt capture < symbols in major element columns
    if contains(colname{c},'_ppm')
        colname{c} = colname{c}(1:end-4);
    elseif contains(colname{c},'_ppb')
        num(i,:) = num(i,:)*1e3;
        colname{c} = colname{c}(1:end-4);
    elseif contains(colname{c},'_ppt')
        num(i,:) = num(i,:)*1e-6;
        colname{c} = colname{c}(1:end-4);
    elseif contains(colname{c},'_per') || contains(colname{c},'_wtper')
        num(i,:) = num(i,:)*1e4;
        if contains(colname{c},'_per')
            colname{c} = colname{c}(1:end-4);
        else
            colname{c} = colname{c}(1:end-6);
        end
    end

    if sum(strcmpi(colname{c},el)) > 0
        colname{c} = [colname{c},'_ppm'];
        flag = 1;
    elseif sum(strcmpi(colname{c},numfield)) > 0 ...
        || sum(strcmpi(colname{c},timefield)) > 0
        flag = 1;
    end

    if flag
        
        lind = contains({txt{i,:}},'bdl') | contains({txt{i,:}},'b.d.l.') ...
            | contains({txt{i,:}},'b.d.l') ...
            | contains({txt{i,:}},'nd') | contains({txt{i,:}},'n.d.') ...
            | contains({txt{i,:}},'bd') | contains({txt{i,:}},'b.d.') ...
            | contains({txt{i,:}},'<dl') | contains({txt{i,:}},'<d.l.');
        num(i,lind) = 0;
        txt(i,lind) = {''};

        % if there is a single '-' remove it, note this does not always 
        % work with en- or em-dash.
        lind = strcmpi({txt{i,:}},'-');
        num(i,lind) = nan;
        txt(i,lind) = {''};
        
        colname{c}
        lind = isempty(num(i,:)) & (contains({txt{i,:}},'<') | contains({txt{i,:}},'>'))
        if any(lind)
            txt{i,lind} = replace({txt{i,lind}},'< ','-');
            txt{i,lind} = replace({txt{i,lind}},'<','-');
            txt{i,lind} = replace({txt{i,lind}},'> ','');
            txt{i,lind} = replace({txt{i,lind}},'>','');
            num(i,lind) = str2num({txt{i,lind}});
        end

%         for j = 1:length(num(i,:))
% 
%             
% 
%             if isnan(num(i,j)) && ~isempty(txt{i,j})
%                 ind = strfind(txt{i,j},'< ');
%                 if ~isempty(ind)
%                     txt{i,j} = ['-',txt{i,j}(ind+1:end)];
% 
%                     [N,ok] = str2num(txt{i,j});
%                     if ok
%                         num(i,j) = N;
%                         continue;
%                     end
%                 end
% 
%                 % subtle difference to above '< ', which includes a space
%                 %Consider using strtrim and regexprep?
%                 ind = findstr(txt{i,j},'<');
%                 if ~isempty(ind)
%                     txt{i,j} = ['-',txt{i,j}(ind+1:end)];
% 
%                     [N,ok] = str2num(txt{i,j});
%                     if ok
%                         num(i,j) = N;
%                         continue;
%                     end
%                 end
% 
%                 % there seem to be so few '>' that we can safely ignore them
%                 ind = findstr(txt{i,j},'>');
%                 if ~isempty(ind)
%                     txt{i,j}(ind) = '';
%                     [N,ok] = str2num(txt{i,j});
%                     if ok
%                         num(i,j) = N;
%                         continue;
%                     end
%                 end
% 
%                 
%             end
%         end
%         % change < to - , exact '-' to NaN, '>' to blank, 'b.d.l.' and 'n.d.' to 0
% 
%         %if sum(strcmpi(colname{c},timefield)) > 0
%         %    % create text-based timescale names
%         %    colname{c} = 
%         %end

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
        for j = 1:nr
            tmp{j,nc+x} = '';
        end
        x = x + 1;
    end
end

data = cell2table(tmp,'VariableNames',colname);

return
