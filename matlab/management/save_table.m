function save_table(data,table_name)
% Outputs MATLAB table format variable to csv file and replaces all NaN in
% numerical columns with blanks to reduce memory
% Inputs: data - table format variable
%         table_name - name of csv file to be saved to

% Write to file
fprintf('Writing %s...\n',table_name)
csv_name = [table_name,'.csv'];
writetable(data,csv_name,'Delimiter',',','QuoteStrings',true);
fid = fopen(csv_name,'rt') ;
X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;
Y = strrep(X, 'NaN', '') ;
clear X
fid2 = fopen(csv_name,'wt') ;
fwrite(fid2,Y) ;
fclose (fid2) ;

end