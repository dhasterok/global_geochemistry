function alkalicratio(data,age_div)
%Ratio of alkalic to subalkalic

for i = 1:length(age_div)-1
    %Alkalic basalt
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & strcmpi(data.rock_type,'alkalic basalt'));
    ab.n(i) = length(ind);
    ab.ind{i} = ind;
    ab.age{i} = data.avg_age(ind);
    
    %Subalkalic basalt
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & strcmpi(data.rock_type,'subalkalic basalt'));
    sb.n(i) = length(ind);
    sb.ind{i} = ind;
    sb.age{i} = data.avg_age(ind);
    
    
    
    %Alkalic gabbro
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & strcmpi(data.rock_type,'alkalic gabbro'));
    ag.n(i) = length(ind);
    ag.ind{i} = ind;
    ag.age{i} = data.avg_age(ind);
    
    %Subalkalic gabbro
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & strcmpi(data.rock_type,'subalkalic gabbro'));
    sg.n(i) = length(ind);
    sg.ind{i} = ind;
    sg.age{i} = data.avg_age(ind);
end

%Just need an age vector of each
alkalicb = data.avg_age(strcmpi(data.rock_type,'alkalic basalt') & data.avg_age >= 0);
subalkalicb = data.avg_age(strcmpi(data.rock_type,'subalkalic basalt') & data.avg_age >= 0);

alkalicg = data.avg_age(strcmpi(data.rock_type,'alkalic gabbro') & data.avg_age >= 0);
subalkalicg = data.avg_age(strcmpi(data.rock_type,'subalkalic gabbro') & data.avg_age >= 0);


figure()
subplot(1,3,1)
histogram(alkalicb,age_div,'DisplayStyle','stairs','EdgeColor','b');
hold on
histogram(subalkalicb,age_div,'DisplayStyle','stairs','EdgeColor','r');
hold off
legend('Alkalic basalt','Subalkalic basalt')
xlabel('Age [Ma]')
ylabel('No. data')
golden

subplot(1,3,2)
histogram(alkalicg,age_div,'DisplayStyle','stairs','EdgeColor','b');
hold on
histogram(subalkalicg,age_div,'DisplayStyle','stairs','EdgeColor','r');
hold off
legend('Alkalic gabbro','Subalkalic gabbro')
xlabel('Age [Ma]')
ylabel('No. data')
golden

subplot(1,3,3)
combinedalkalic = vertcat(alkalicb,alkalicg);
combinedsubalkalic = vertcat(subalkalicb,subalkalicg);
histogram(combinedalkalic,age_div,'DisplayStyle','stairs','EdgeColor','b');
hold on
histogram(combinedsubalkalic,age_div,'DisplayStyle','stairs','EdgeColor','r');
hold off
legend('Alkalic (combined)','Subalkalic (combined)')
xlabel('Age [Ma]')
ylabel('No. data')
golden


figure()
subplot(1,2,1)
ratiob = (ab.n)./(sb.n);
h = bar(age_div(1:end-1),(ab.n + ag.n)./((ab.n + ag.n)+(sb.n + sg.n)),'histc');
set(h,'FaceColor',[0.5 0.5 0.5])

xlabel('Age [Ma]')
ylabel('Alkalic ratio')
golden

subplot(1,2,2)
h = bar(age_div(1:end-1),(sb.n + sg.n)./((ab.n + ag.n)+(sb.n + sg.n)),'histc');
set(h,'FaceColor',[0.5 0.5 0.5])

xlabel('Age [Ma]')
ylabel('Subalkalic ratio')
golden


[age_div(1:end-1)' ab.n' sb.n' ag.n' sg.n']


figure()
combineda = (ag.n+ab.n);
combineds = (sg.n+sb.n);
combinedratio = combineda./(combineda + combineds);

x = (age_div(2:end) + age_div(1:end-1))./2;
x = vertcat(ones(1,size(x,2)),x);

size(combinedratio')
size(x')


eqval = regress(combinedratio',x');
slope = eqval(2);
intercept = eqval(1);

%Shift them all to 1500 Ma equivalent slope - x(8)

for i = 1:length(x)
    combinedratio(1,i) = combinedratio(1,i).*(x(i)-x(8)).*slope;
end



h = bar(age_div(1:end-1),combinedratio(1,:))


return