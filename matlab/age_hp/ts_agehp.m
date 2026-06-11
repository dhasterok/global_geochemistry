function output = ts_agehp(data,division)
%Division: eon,era,period,epoch,age


data = data(:,{'sample_id','rock_name','time_period','age','heat_production','country','rock_type'});
data = data(data.heat_production >= 0,:);
data = data(~strcmpi(data.time_period,'') | ~isnan(data.age),:);

ts = readtable('ts.csv');
%Make sure to check for:
%Different terms: upper/late, middle/mid,lower/early
%Dashed: upper-/late-, middle-/mid-, lower-/early-
%Different names: e.g. upper/middle/lower permian instead of actual names
%Hyphenated ranges: e.g. early ordovician-late cambrian and variants
%Different spellings: Archean/Archaean, paleozoic/palaeozoic

%Whatever resolution picked, make every subdivision equal to the resolution
%wanted e.g. upper-triassic -> triassic

%Check that division is correct

division = lower(division);
if ~ismember(division, {'eon','era','period','epoch','age'})
    error('Incorrect division term');
    return;
end

%1. Replace all '-' with ' '
data.time_period = regexprep(data.time_period,'-',' ');
data.time_period = regexprep(data.time_period,' ','_');
ts.name = regexprep(ts.name,' ','_');

%2. Replace all archean/paleozoic etc with australian version
data.time_period = regexprep(data.time_period,'archean','archaean','ignorecase');
data.time_period = regexprep(data.time_period,'paleo','palaeo','ignorecase');

%3. Compare the string to known intervals from ts
%Replace all late/mid/early with upper/lower/middle
%Can include more checks, but for now this will do
data.time_period = regexprep(data.time_period,'late','upper','ignorecase');
%data.time_period = regexprep(data.time_period,'mid','middle','ignorecase');
data.time_period = regexprep(data.time_period,'early','lower','ignorecase');

%3.1 Replace all incorrectly named time periods that i know of
data.time_period = regexprep(data.time_period,'middledle','middle','ignorecase');
data.time_period = regexprep(data.time_period,'middle cambrian','cambrian','ignorecase');
data.time_period = regexprep(data.time_period,'siluro devonian','palaeozoic','ignorecase');
data.time_period = regexprep(data.time_period,'cambro ordovician','palaeozoic','ignorecase');
data.time_period = regexprep(data.time_period,'lower precambrian','archaean','ignorecase');
data.time_period = regexprep(data.time_period,'palaeo mesoarchaean','archaean','ignorecase');
data.time_period = regexprep(data.time_period,'upper neoarchaean','neoarchaean','ignorecase');
data.time_period = regexprep(data.time_period,'middle neoproterozoic','neoproterozoic','ignorecase');
data.time_period = regexprep(data.time_period,'upper neoproterozoic','neoproterozoic','ignorecase');
data.time_period = regexprep(data.time_period,'upper tertiary','miocene','ignorecase');




%5. Find indices where time_period is in the ts file
h = waitbar(0,'Initializing waitbar...');

ts_ind = [];

for i = 1:size(data,1)
    if ~isnan(data.age(i))
        ts_ind = [ts_ind i];
        perc = i/size(data,1);
        waitbar(perc,h,sprintf('Searching for time periods \n %0.2f%%',perc*100))
        continue;
    elseif ~strcmpi(data.time_period(i),'')
        ts_ind = [ts_ind i];
        perc = i/size(data,1);
        waitbar(perc,h,sprintf('Searching for time periods \n %0.2f%%',perc*100))
        continue;
    end
end
close(h)

%5. Division fixes for these indices
%Use age ranges

epoch_ind = find(strcmpi(ts.type,'epoch'));
period_ind = find(strcmpi(ts.type,'period'));
era_ind = find(strcmpi(ts.type,'era'));
eon_ind = find(strcmpi(ts.type,'eon'));

h = waitbar(0,'Initializing waitbar...');
if strcmpi(division,'age')
    %All era/period/epoch/age within an eon get rounded up
    
elseif strcmpi(division,'epoch')
    %All period/epoch/age within an era get rounded up, eon ignored
    

    
elseif strcmpi(division,'period')
    %All epoch/age within a period get rounded up, eon/era ignored
    
    
elseif strcmpi(division,'era')
    %All age within an epoch get rounded up, eon/era/period ignored
    for i = 1:length(ts_ind)
        %If is not an age, delete the reference
        if ~isnan(data.age(ts_ind(i)))
            continue;
        end
        if ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(strcmpi(ts.type,'age')))) || ...
                ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(strcmpi(ts.type,'epoch')))) || ...
                ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(strcmpi(ts.type,'period'))))
            for j = 1:length(era_ind)
                if ts.age_lower(era_ind(j)) <= ts.age_lower(strcmpi(data.time_period(ts_ind(i)),ts.name)) ...
                        && ts.age_upper(era_ind(j)) >= ts.age_upper(strcmpi(data.time_period(ts_ind(i)),ts.name))
                    data.time_period(ts_ind(i)) = ts.name(era_ind(j));
                    
                    break;
                end
            end
        elseif ~ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(find(strcmpi(ts.type,'era')))))
            data.time_period(ts_ind(i)) = {''};
        end
        perc = i/length(ts_ind);
        waitbar(perc,h,sprintf('Setting to correct age divisions \n %0.2f%%',perc*100))
    end
    close(h)
    
elseif strcmpi(division,'eon')
    %Eon/era/period/epoch ignored
    for i = 1:length(ts_ind)
        %If is not an age, delete the reference
        if ~isnan(data.age(ts_ind(i)))
            continue;
        end
        if ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(strcmpi(ts.type,'age')))) || ...
                ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(strcmpi(ts.type,'epoch')))) || ...
                ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(strcmpi(ts.type,'period')))) || ...
                ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(strcmpi(ts.type,'era'))))
            for j = 1:length(eon_ind)
                if ts.age_lower(eon_ind(j)) <= ts.age_lower(strcmpi(data.time_period(ts_ind(i)),ts.name)) ...
                        && ts.age_upper(eon_ind(j)) >= ts.age_upper(strcmpi(data.time_period(ts_ind(i)),ts.name))
                    data.time_period(ts_ind(i)) = ts.name(eon_ind(j));
                    
                    break;
                end
            end
        elseif ~ismember(lower(data.time_period(ts_ind(i))),lower(ts.name(find(strcmpi(ts.type,'eon')))))
            data.time_period(ts_ind(i)) = {''};
        end
        perc = i/length(ts_ind);
        waitbar(perc,h,sprintf('Setting to correct age divisions \n %0.2f%%',perc*100))
    end
    close(h)
    
end

ts_ind = ts_ind(~strcmpi(data.time_period(ts_ind),'') | ~isnan(data.age(ts_ind)));

size(data(isnan(data.age(ts_ind,:)) & strcmpi(data.time_period(ts_ind),'neoproterozoic'),:),1)

%__________________________________________________________________________

%         Calculate the average/median HP as boxplots in y only
%__________________________________________________________________________

%Based on eon/epoch/period pick up top, the age ranges to use
%Search through time_period for the name if no age, and if age use that
%Combine the indices: these combined are the bin to feed to statistics

%Data is pre-filtered and ONLY epoch/period etc specified at start now

%Loop through iterative bplot call with set width. X position is centre of
%box, width is equal to total width of box wanted
%bplot([data],x position,'width',width_value)
if strcmpi(division,'age')
    tstable_ind = find(strcmpi(ts.type,'age'));
elseif strcmpi(division,'epoch')
    tstable_ind = find(strcmpi(ts.type,'epoch'));
elseif strcmpi(division,'period')
    tstable_ind = find(strcmpi(ts.type,'period'));
elseif strcmpi(division,'era')
    tstable_ind = find(strcmpi(ts.type,'era'));
elseif strcmpi(division,'eon')
    tstable_ind = find(strcmpi(ts.type,'eon'));
end

%Create an empty table for each time period
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);

h = waitbar(0,'Initializing waitbar...');

for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
        
        
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)


%__________________________________________________________________________

%          Plot the time period box + the actual age dates
%__________________________________________________________________________


%After compare, if isnt found BUT first part of is upper/middle/lower and
%second part is found within - Different names: e.g. upper/middle/lower 
%permian instead of actual names

%Loop through iterative bplot call with set width. X position is centre of
%box, width is equal to total width of box wanted
%bplot([data],x position,'width',width_value)

%REWRITE THIS to save ind into the array rather than HP values so we can
%extract info on them if needed. Will also allow removal of all these
%repeating loops.

%ADD NUMBER OF DATA IN EACH BIN

figure()
subplot(1,2,1)
title('All data')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

%Create an empty table for each time period
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);

h = waitbar(0,'Initializing waitbar...');

for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
        
        
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)

subplot(1,2,2)
title('Without ALL Aus')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

output = databins;


figure()
%Granite/rhyolite
subplot(3,2,1)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'granite') || strcmpi(data.rock_type(ts_ind(i)),'rhyolite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'granite') || strcmpi(data.rock_type(ts_ind(i)),'rhyolite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Granite/Rhyolite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')


%Dacite/Granodiorite
subplot(3,2,2)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'dacite') || strcmpi(data.rock_type(ts_ind(i)),'granodiorite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'dacite') || strcmpi(data.rock_type(ts_ind(i)),'granodiorite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Dacite/Granodiorite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

%Diorite/Andesite
subplot(3,2,3)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'diorite') || strcmpi(data.rock_type(ts_ind(i)),'andesite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'diorite') || strcmpi(data.rock_type(ts_ind(i)),'andesite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Diorite/andesite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

%Basaltic andesite/gabbroic diorite
subplot(3,2,4)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'basaltic andesite') || strcmpi(data.rock_type(ts_ind(i)),'gabbroic diorite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'basaltic andesite') || strcmpi(data.rock_type(ts_ind(i)),'gabbroic diorite'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Basaltic andesite/gabbroic diorite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

%Subalkalic basalt/subalkalic gabbro
subplot(3,2,5)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'subalkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'subalkalic gabbro'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'subalkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'subalkalic gabbro'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Subalkalic basalt/subalkalic gabbro')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')


%Alkalic basalt/Alkalic gabbro
subplot(3,2,6)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'alkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'alkalic gabbro'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'alkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'alkalic gabbro'))
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('alkalic basalt/alkalic gabbro')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')


%-------------
%Without Aus
%-------------
figure()
%Granite/rhyolite
subplot(3,2,1)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'granite') || strcmpi(data.rock_type(ts_ind(i)),'rhyolite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'granite') || strcmpi(data.rock_type(ts_ind(i)),'rhyolite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Granite/Rhyolite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')


%Dacite/Granodiorite
subplot(3,2,2)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'dacite') || strcmpi(data.rock_type(ts_ind(i)),'granodiorite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'dacite') || strcmpi(data.rock_type(ts_ind(i)),'granodiorite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Dacite/Granodiorite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

%Diorite/Andesite
subplot(3,2,3)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'diorite') || strcmpi(data.rock_type(ts_ind(i)),'andesite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'diorite') || strcmpi(data.rock_type(ts_ind(i)),'andesite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Diorite/andesite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

%Basaltic andesite/gabbroic diorite
subplot(3,2,4)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'basaltic andesite') || strcmpi(data.rock_type(ts_ind(i)),'gabbroic diorite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'basaltic andesite') || strcmpi(data.rock_type(ts_ind(i)),'gabbroic diorite')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Basaltic andesite/gabbroic diorite')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')

%Subalkalic basalt/subalkalic gabbro
subplot(3,2,5)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'subalkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'subalkalic gabbro')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'subalkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'subalkalic gabbro')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('Subalkalic basalt/subalkalic gabbro')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')


%Alkalic basalt/Alkalic gabbro
subplot(3,2,6)
databins  = cell(2,length(tstable_ind));
databins(1,:) = ts.name(tstable_ind,:);
h = waitbar(0,'Initializing waitbar...');
for i = 1:length(ts_ind)
    if ~isnan(data.age(ts_ind(i)))
        for j = 1:length(tstable_ind)
            if data.age(ts_ind(i)) >= ts.age_lower(tstable_ind(j)) && ...
                    data.age(ts_ind(i)) <= ts.age_upper(tstable_ind(j)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'alkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'alkalic gabbro')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    elseif isnan(data.age(ts_ind(i))) && ~strcmp(data.time_period(ts_ind(i)),'')
        for j = 1:length(tstable_ind)
            if strcmpi(databins(1,j),data.time_period(ts_ind(i),:)) && ...
                    (strcmpi(data.rock_type(ts_ind(i)),'alkalic basalt') || strcmpi(data.rock_type(ts_ind(i)),'alkalic gabbro')) && ...
                    ~strcmpi(data.country(ts_ind(i)),'australia') && ...
                    ~strcmpi(data.country(ts_ind(i)),'aus') && ...
                    ~strcmpi(data.country(ts_ind(i)),'au')
                databins{2,j} = [databins{2,j} data.heat_production(ts_ind(i))];
                break;
            end
        end
    end
    perc = i/length(ts_ind);
    waitbar(perc,h,sprintf('Partitioning up HP values for timescales\n %0.2f%%',perc*100))
end
close(h)
title('alkalic basalt/alkalic gabbro')
hold on
for i = 1:size(databins,2)
    bplot(databins{2,i},((ts.age_lower(strcmpi(ts.name,databins{1,i})))+...
        (ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,...
        'width',(ts.age_upper(strcmpi(ts.name,databins{1,i})))-(ts.age_lower(strcmpi(ts.name,databins{1,i}))));
    text(((ts.age_lower(strcmpi(ts.name,databins{1,i})))+(ts.age_upper(strcmpi(ts.name,databins{1,i}))))/2,3,int2str(length(databins{2,i})))
end
hold off
xlabel('Age [Mya]')
ylabel('Heat production [\muW m^{-3}]')


return