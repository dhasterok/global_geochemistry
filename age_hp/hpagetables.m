max_id = max(data.spreadsheet_id);
spreadsheet_table = cell(max_id,2);
spreadsheet_table(:,2) = {0};
spreadsheet_table(:,1) = num2cell([1:max_id]');

h = waitbar(0,'Initializing waitbar...');
for i = 1:size(data,1)
    if (i./size(data,1))*100 > 99.60
        i
    end
    spreadsheet_table{data.spreadsheet_id(i),2} = spreadsheet_table{data.spreadsheet_id(i),2}+1;
    perc = 100*(i/size(data,1));
    waitbar(perc/100,h,sprintf('Progress:\n %0.2f%%',perc))
end
close(h)
sortrows(spreadsheet_table,-2)