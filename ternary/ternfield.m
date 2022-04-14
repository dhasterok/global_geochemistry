function [A1,A2,A3] = ternfield(data,verticies,wt2mol,varargin)
% [A1,A2,A3] = ternfield(data,verticies,wt2mol,varargin)

% Errorlog:
%   26 Feb 2019 - fixed an error in feo_tot to fe2o3_tot conversion

a1 = verticies{1};
a2 = verticies{2};
a3 = verticies{3};

if any(strcmp('Fe2O3',verticies))    
    cf = molecularwt('Fe2O3')/(2*molecularwt('FeO'));
    data.fe2o3_tot = data.feo_tot*cf;
end

% check number of inputs and create an initial index of logicals if
% necessary
field = 'heat_production';
dflag = 0;
rflag = 0;
pflag = 0;
indp = logical(ones([height(data) 1]));
if nargin > 3
    c = 1;
    while c <= nargin - 3
        switch lower(varargin{c})
            case 'field'
                field = varargin{c+1};
                if ~any(strcmpi(field,data.Properties.VariableNames))
                    error('Field not found in input data.');
                end
                c = c + 2;
            case 'difference'
                dflag = 1;
                indp = varargin{c+1};
                indq = varargin{c+2};
                c = c + 3;
            case 'rtype'
                rflag = 1;
                rfield = varargin{c+1};
                rtypes = varargin{c+2};
                c = c + 3;
            case 'polyplot'
                pflag = 1;
                ptype = varargin{c+1};
                c = c + 2;
            otherwise
                error('Unknown option.');
        end
        
    end
end

% select only values greater than 0
ind0 = data{:,field} > 0;

% if field is heat production, take log of data
if strcmp(field,'heat_production')
    data.heat_production = log10(data.heat_production);
end

% parse data
A1 = parsea(data,a1,wt2mol);
A2 = parsea(data,a2,wt2mol);
A3 = parsea(data,a3,wt2mol);

% sum of data
T = sum(A1 + A2 + A3,2);

% normalize data to total
ind = ~isnan(T) & ind0;
A1 = A1(ind)./T(ind);
A2 = A2(ind)./T(ind);
A3 = A3(ind)./T(ind);

data = data(ind,:);
indp = indp(ind);
if dflag
    indq = indq(ind);
end

%figure;
%if sum(A1 < 0) > 0
%    ternextend(0.5);
%end
%hold on;
%ternary(a1,a2,a3);
%t = ternplot(A1,A2,A3,log10(data.heat_production(ind)));
%set(t,'SizeData',2*ones(size(get(t,'SizeData'))));
%caxis([-1 1]);
%cbar('log_{10} A');

figure;
ia = logical(ones([height(data) 1]));
if dflag
    ip = ia & indp;
    iq = ia & indq;
    ia = ip | iq;
end

if sum(A1 < 0) > 0
    if dflag
        outp = ternsurf(A1(ip),A2(ip),A3(ip),data{ip,field},0.05,0.5);
        outq = ternsurf(A1(iq),A2(iq),A3(iq),data{iq,field},0.05,0.5);
    end
    outa = ternsurf(A1(ia),A2(ia),A3(ia),data{ia,field},0.05,0.5);
else
    if dflag
        outp = ternsurf(A1(ip),A2(ip),A3(ip),data{ip,field},0.05);
        outq = ternsurf(A1(iq),A2(iq),A3(iq),data{iq,field},0.05);
    end
    outa = ternsurf(A1(ia),A2(ia),A3(ia),data{ia,field},0.05);
end
if dflag & rflag
    subplot(6,3,[4,7]);
elseif dflag & ~rflag
    subplot(2,3,1);
else
    subplot(2,2,1);
end
trisurf(outa.tri,outa.xv,outa.yv,outa.median);
hold on;
%shading interp;
shading flat;
set(gca,'View',[0 90]);

ternary(a1,a2,a3);

if strcmp(field,'heat_production')
    caxis([-1 1]);
    cbar('log_{10} A');
else
    cbar('median');
end
%colormap(hot);

if dflag & rflag
    subplot(6,3,[10,13]);
elseif dflag & ~rflag
    subplot(2,3,4);
else
    subplot(2,2,3);
end
trisurf(outa.tri,outa.xv,outa.yv,outa.sd);
%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
ternary(a1,a2,a3);
caxis([0 2]);
cbar('\sigma_{\mu}');

if dflag & rflag
    if rflag
        subplot(6,3,[6,9]);
    else
        subplot(2,3,2);
    end
    trisurf(outp.tri,outp.xv,outp.yv,outp.nbin);
    cbar('No. data1');
else
    subplot(2,2,2);
    trisurf(outa.tri,outa.xv,outa.yv,outa.nbin);
    cbar('No. data');
end

%shading interp;
shading flat;
set(gca,'View',[0 90]);
hold on;
ternary(a1,a2,a3);
caxis([0 3]);

if dflag
    if rflag
        subplot(6,3,[12,15]);
    else
        subplot(2,3,5);
    end    
    trisurf(outq.tri,outq.xv,outq.yv,outq.nbin);
    %shading interp;
    shading flat;
    set(gca,'View',[0 90]);
    hold on;
    ternary(a1,a2,a3);
    caxis([0 3]);
    cbar('No. data2');
    %cbar('\sigma_{\mu}(N-1)^{-1/2}');

    if rflag
        subplot(6,3,[14,17]);
    else
        subplot(2,3,3);
    end     
    trisurf(outp.tri,outp.xv,outp.yv, outp.median - outq.median);
    %shading interp;
    shading flat;
    set(gca,'View',[0 90]);
    hold on;
    ternary;
    ternary(a1,a2,a3);
    caxis([-1 1]);
    %colormap(rwb);
    cbar('log_{10} sed. - log_{10} metased.');
    
    if rflag
        subplot(6,3,[8,11]);
    else
        subplot(2,3,6);
    end    
    histogram(outp.median - outq.median,'BinEdges',[-2:0.2:2]);
    xlim([-2 2]);
    xlabel('P - Q');
    golden;
end

if rflag
    if dflag
        subplot(6,3,[2,5]);
    else
        subplot(2,2,4);
    end
    ternary(a1,a2,a3);
    c = 1;
    %map = colormap(jet(length(rtypes)));
    for i = 1:length(rtypes)
        ind = strcmp(rtypes{i},data{:,rfield}) & ia;
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
    for i = 1:length(A1m)
        u(i) = ternplot(A1m(i),A2m(i),A3m(i),'^');
    end
    set(u,'MarkerFaceColor',[1 1 1]);

    % text for rock types
    %[x,y] = tern2xy(A1m,A2m,A3m);
    %for i = 1:length(x)
    %    text(x(i),y(i),rtypes{tind(i)});
    %end
    % legend for rock types
    legend(u,rtypes{tind});
    
    if pflag
        if strcmp(ptype,'sed');
            % load QFL polygons
            sedgons = load_sedgons;
            nsed = length(sedgons(:,1));

            for i = 1:nsed
                h = ternplot(sedgons{i,1}(:,1),sedgons{i,1}(:,2),sedgons{i,1}(:,3),'k-');
                set(h,'LineWidth',0.25);
            end
        end
        %elseif strcmp(ptype,'ig');
    end
end

return