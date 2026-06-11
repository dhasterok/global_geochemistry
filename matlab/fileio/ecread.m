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
            fields{i} = [fields{i},'_',str];
        else
            fields{i} = str;
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
            init = struct(fields{i},[]);
        else
            init = setfield(init,fields{i},[]);
        end
    else                % text
        if i == 1
            init = struct(fields{i},'');
        else
            init = setfield(init,fields{i},'');
        end
    end
end

for j = 2:nrow
    % copy empty initialized structure to next data array
    temp = init;
    % fill data array with values
    for i = 1:ncol
        if numfield(i)
            temp = setfield(temp,fields{i},alldata{j,i});
        else
            if isnumeric(alldata{j,i})
                temp = setfield(temp,fields{i},num2str(alldata{j,i}));
            else
                temp = setfield(temp,fields{i},alldata{j,i});
            end
        end
    end
    data(j-1) = temp;
end

% check K vs K2O
% if ~isfield(data(i),'K2O') & ~isfield(data(i),'K')
%     return;
% end
% for i = 1:length(data)
%     if isempty(data(i).K)
%         continue;
%     end
%     K = data(i).K;
%     if ~isempty(data(i).K2O)
%         Ko = data(i).K2O * 2*39.0986/(2*39.0986 + 15.9994)*10000;
%         if K/Ko < 0.1
%             K = K*10000;
%             if K/Ko > 10
%                 warning(['Oops problem with K in entry ',i]);
%                 continue;
%             else
%                 data(i).K = K;
%             end
%         end
%     elseif K < 100
%         data(i).K = K*10000;
%     end    
% end

return
