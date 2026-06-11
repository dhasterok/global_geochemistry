function a = sqwavefill(whisker,age,age_div,varargin)
%whisker is an nx5 matrix
%5,25,50,75,95 percentile y positions
%age matrix (size n), the corresponding x position of the whiskers
%age_div used
%varargin -> color as string e.g. 'b'


lx = length(age);
ly = size(whisker,1);

X = zeros(2*lx,1);
Y = zeros(2*ly,5);

X(1:2:2*lx-1) = age_div(1:lx);
X(2:2:2*lx-2) = age_div(2:lx);
X(end)=age_div(end);

Y(1:2:2*ly-1,:) = whisker(1:ly,:);
Y(2:2:2*ly,:)   = whisker(1:ly,:);
    
if nargin>3
    color = varargin{1};
else
    color = 'b';
end
if nargin>4
    dname = varargin{2};
else
    dname = 'fill';
end

hold on
%5 95
limits_plot = plot(X,Y(:,[1 5]),'Color',color);


%Quartile fill
ind = find(isnan(Y(:,2)) | isnan(Y(:,4)));
if isempty(ind)
    a = fill([X;flipud(X)],([Y(:,2);flipud(Y(:,4))]),color,'DisplayName',dname,'FaceAlpha',0.3,'EdgeColor','None');
else
    if ind(1) > 1
        ind = [0; ind];
    end
    if ind(end) < size(Y,1)
        ind = [ind; size(Y,1)+1];
    end
    
    for i = 1:length(ind)-1
        ii = [ind(i)+1:ind(i+1)-1];
        x = X(ii);
        y = Y(ii,[2 4]);
        a = fill([x;flipud(x)],([y(:,1);flipud(y(:,2))]),color,'DisplayName',dname,'FaceAlpha',0.3,'EdgeColor','None');
    end
end

%Median
plot(age,whisker(:,3),'Color',color,'Marker','.','MarkerSize',8,'linestyle', 'none');

median_plot = plot(X,Y(:,3),'Color',color,'linewidth',0.5);

hold off
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');

return