function A = parsea(data,atxt,wt2mol)

% find locations of operators
op = sort([strfind(atxt,'*'),strfind(atxt,'/'), ...
    strfind(atxt,'+'),strfind(atxt,'-')]);
Nop = length(op);
% if no operators, get data, convert(?) to moles, and return
if Nop == 0    
    if sum(strcmpi(atxt,data.Properties.VariableNames)) == 1 | ...
                strcmpi(atxt,'FeO') | strcmpi(atxt,'Fe2O3') | ...
                strcmpi(atxt,'CaO_adj') 
        % get data from field
        if strcmpi(atxt,'FeO')
            A = data{:,'feo_tot'};
        elseif strcmpi(atxt,'Fe2O3')
            A = data{:,'fe2o3_tot'};
        elseif strcmpi(atxt,'CaO_adj')
            A = data{:,'cao'} ...
                - 10/3*data{:,'p2o5'}*molecularwt('CaO')/molecularwt('P2O5');
            A(A<0) = 0;
            atxt = 'CaO';
        else
            %A = table2array(data(:,lower(atxt)));
            if any(strcmp(atxt,data.Properties.VariableNames))
                A = data{:,atxt};
            else
                A = data{:,lower(atxt)};
            end
        end        
        % remove values less than zero
        A(A<0) = NaN;
        
        % convert to moles if necessary
        if wt2mol == 1
            A = A/molecularwt(atxt);
        end
    else
        error(['Field ',atxt,' not found.']);
    end
    return 
end

% grab all fields first and place into array by columns.
istart = 1;
for i = 1:Nop+1
    if i == Nop+1;
        iend = length(atxt);
    else
        iend = op(i)-1;
    end
    
    str = strtrim(atxt(istart:iend));
    [num,ok] = str2num(str);
    if ok
        x(:,i) = num*ones([length(data.sio2),1]);
    else
        if any(strcmpi(str,data.Properties.VariableNames)) | ...
                strcmpi(str,'FeO') | strcmpi(str,'Fe2O3') | ...
                strcmpi(str,'CaO_adj')
            % get data from field
            if strcmpi(str,'FeO')
                x(:,i) = data{:,'feo_tot'};
            elseif strcmpi(str,'Fe2O3')
                x(:,i) = data{:,'fe2o3_tot'};
            elseif strcmpi(str,'CaO_adj')
                x(:,i) = data{:,'cao'} ...
                    - 10/3*data{:,'p2o5'}*molecularwt('CaO')/molecularwt('P2O5');
                x(x(:,i)<0) = 0;
                str = 'CaO';
            else
                if any(strcmp(str,data.Properties.VariableNames))
                    x(:,i) = data{:,str};
                else
                    x(:,i) = data{:,lower(str)};
                end
                %x(:,i) = table2array(data(:,str));
            end
            % remove zero values from field
            x(x(:,i)<0) = NaN;
            
            % convert to moles if necessary
            if wt2mol == 1
                x(:,i) = x(:,i)/molecularwt(str);
            end
            
        else
            error(['Field ',str,' not field.']);
        end
    end
    istart = iend + 2;
end


Nx = Nop + 1;
% Order of operations / first, then *, then + and - last.
% Convert divisions to multiplications first by taking reciporical of
% vecton first and change / to * in operator list
for i = 1:Nop
    if strcmp(atxt(op(i)),'/')
        x(:,i+1) = 1./x(:,i+1);
        atxt(op(i)) = '*';
    end
end

% Perform all * and / operations.  For each operation, 
% (1) store in lower position
% (2) move shift all columns down one
% (3) remove operator from list

for i = Nop:-1:1
    m = atxt(op(i));
    if strcmp(m,'*')
        x(:,i) = colmath(x(:,i),x(:,i+1),m);
        
        % shift all columns above i + 2 down by 1
        if i + 2 < Nx
            x = x(:,[1:i,i+2:Nx]);
        end
        Nx = Nx - 1;
        
        % delete operator index
        if i ~= Nop
            op = op([1:i-1 i+1:end]);
        else
            op = op([1:i-1]);
        end
        Nop = Nop - 1;
    end
end

% Operations + and - are simpler.  Always replace first column.
for i = 1:Nop
    m = atxt(op(i));
    x(:,1) = colmath(x(:,1),x(:,i+1),m);
end

% Final result should be 1st column
A = x(:,1);

return

function x = colmath(x1,x2,op)

switch op
    case '*'
        x = x1 .* x2;
    case '/'
        x = x1 ./ x2;
    case '+'
        x = x1 + x2;
    case '-'
        x = x1 - x2;
end

return