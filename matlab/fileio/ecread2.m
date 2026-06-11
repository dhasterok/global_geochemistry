function data = ecread(filename);
% ECREAD - Reads Earthchem xls output
%
%   DATA = ECREAD(FILENAME) where FILENAME is an excel spreadsheet containing
%   earthchem.org data.  The output, DATA, is a structure array containing the
%   earthchem data with each field being a column from the original.

[~,~,alldata] = xlsread(filename);

[nrow,ncol] = size(alldata);
% field types: (0) text; (1) number
numfield = logical(ones([1,ncol]));
for i = 1:ncol
    % construct field name
    temp = alldata{1,i};
    if ~ischar(temp)
        temp = num2str(temp);
    end
    fields{i} = '';
    flag = 0;
    while ~isempty(temp)
        [str,temp] = strtok(temp);
        if flag == 1
            fields{i} = [fields{i},'_',upper(str)];
        else
            fields{i} = upper(str);
            flag = 1;
        end
    end

    % determine if the field should be text or numerical
    for j = 2:nrow
        if ~isempty(alldata{j,i}) & ischar(alldata{j,i})
            str = alldata{j,i};

            % test for below detection limits
            % just because alldata includes text does not imply a text field.
            % some chemical data have a qualifier '<' preceeding.  this means
            % the data fall below the reported detection limit.  convert '<' to
            % a '-' so that it can be converted to a number.  since negative
            % concentrations cannot exist, filtering for negative numbers can be
            % used to identify these.
            if strcmp(str(1),'<')
                str(1) = '-';
            end

            % ignore '>' preceeding, there generally aren't enough of these in
            % the entire database to skew the results
            if strcmp(str(1),'>')
                str(1) = '';
            end

            % check and convert to number and adjust alldata
            [x,ok] = str2num(str);
            if ok
                alldata{j,i} = x;
                continue;
            end

            numfield(i) = 0;
            break;
        end
    end

    % initialize the field as empty
    if numfield(i)    % numerical
        if i == 1
            data = struct(fields{i},nan(nrow-1,1));
        else
            data = setfield(data,fields{i},nan(nrow-1,1));
        end
    else                % text
        if i == 1
            data = struct(fields{i},{cell(nrow-1,1)});
        else
            data = setfield(data,fields{i},cell(nrow-1,1));
        end
    end
end

% fill data array with values
for i = 1:ncol
    if numfield(i)
        temp = nan(nrow-1,1);
        for j = 2:nrow
            if ~isempty(alldata{j,i})
                temp(j-1,1) = alldata{j,i};
            end
        end
    else
        temp = cell(nrow-1,1);
        for j = 2:nrow
            if isnumeric(alldata{j,i})
                temp{j-1,1} = num2str(alldata{j,i});
            else
                temp{j-1,1} = alldata{j,i};
            end     
        end
    end
    data = setfield(data,fields{i},temp);
    data.len = nrow-1;
    clear temp;
end

return
