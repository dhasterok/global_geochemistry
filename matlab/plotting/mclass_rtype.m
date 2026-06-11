function [A1,A2,A3] = mclass_rtype(data,verticies,wt2mol,field,types)

a1 = verticies{1};
a2 = verticies{2};
a3 = verticies{3};
% check number of inputs and create an initial index of logicals if
% necessary

% parse data
A1 = parsea(data,a1,wt2mol);
A2 = parsea(data,a2,wt2mol);
A3 = parsea(data,a3,wt2mol);
%B = A2 - (data.NA2O/molecularwt('Na2O') + data.K2O/molecularwt('K2O'));
%[min(B) max(B) mean(B(~isnan(B))) std(B(~isnan(B)))]

% Create AFM diagram
T = sum(A1 + A2 + A3,2);

ind = ~isnan(T) & T > 0;
A1 = A1(ind)./T(ind);
A2 = A2(ind)./T(ind);
A3 = A3(ind)./T(ind);

data = data(ind,:);

ternary(a1,a2,a3);
c = 1;
%map = colormap(jet(length(types)));
for i = 1:length(types)
    ind = strcmp(types{i},data{:,field});
    if sum(ind) == 0
        continue;
    end
    
    %if i == 1;
    [x,y] = tern2xy(A1(ind),A2(ind),A3(ind));
    [N,C] = hist3([x,y],[50 50]);
    hold on;
        
    n = N(N~=0);
    n = sort(n(:));
    
    %figure;
    %plot(n);
    n10per = 0.1*sum(n(:));         % 10 percent of data
    itmp = sum(cumsum(n) < n10per); % index of 10 percent
    level = n(itmp);                % contour level encompasing 90% of the data
    %hold on;
    %plot([0 length(n(:))],[level level],'-');
    
    %figure;
    %t = ternplot(A1(ind),A2(ind),A3(ind),'.');
    %set(t,'Color',map(i,:));
    
    %contour(C{1},C{2},N','r-');
    contour(C{1},C{2},N',[level level],'-');
        %figure;
        %hist(N/length(x));
    %end

    A1m(c) = median(A1(ind));
    A2m(c) = median(A2(ind));
    A3m(c) = median(A3(ind));
    
    tind(c) = i;
    
    c = c + 1;
end

if sum(A1m < 0) > 0
    ternextend(0.5);
end
hold on;
t = ternplot(A1m(1),A2m(1),A3m(1),'p');
set(t,'MarkerFaceColor',[1 1 1]);
u = ternplot(A1m(2:end),A2m(2:end),A3m(2:end),'^');
set(u,'MarkerFaceColor',get(t,'MarkerEdgeColor'));

[x,y] = tern2xy(A1m,A2m,A3m);
for i = 1:length(x)
    text(x(i),y(i),types{tind(i)});
end
legend(types{tind});

return