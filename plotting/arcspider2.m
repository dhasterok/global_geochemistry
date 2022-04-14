function [t,ellist,data] = arcspider(arcs,data,varargin)
% ARCSPIDER - Computes spider diagrams
%
%   [arcs,ellist,data] = arcspider(arcs,data,varargin)
%
%   Input options:
%       SortField
%       ElementList
%       NormRef         Reference for normalizing spider diagram, use
%                       'mean' to normalize to the mean of the each
%       FieldLabel
%       RefField
%       RefVal
%       SubSet          Logical index

% elements and reference model
normref = 'G13';
ellist = {'Cs', 'Rb', 'Ba', 'Th', ...
    'U', 'Nb', 'Ta', 'P', 'K', ...
    'La', 'Ce', 'Pb', 'Mo', ...
    'Pr', 'Sr', 'Ga', 'Zr', 'Hf', ...
    'Nd', 'Sm', 'Eu', 'Li', ...
    'Ti', 'Gd', 'Tb', 'Dy', ...
    'Y', 'Ho', 'Er', ...
    'Tm', 'Yb', 'Lu', 'Zn', ...
    'Mn', 'V', 'Sc', 'Co', ...
    'Cu', 'Ni', 'Cr'};

% convert K2O to K
data.k_ppm = data.k2o*2*molecularwt('K')/molecularwt('K2O')*1e4;

sind = true(height(data),1);
reffield = 'sio2';
if nargin > 2
    opt = 1;
    while opt + 2 < nargin
        switch lower(varargin{opt})
            case 'sortfield'
                sortfield = varargin{opt+1};
            case 'fieldlabel'
                fieldlbl = varargin{opt+1};
            case 'reffield'
                reffield = varargin{opt+1};
            case 'refval'
                refval = varargin{opt+1};
            case 'subset'
                sind = varargin{opt+1};
            case 'subtext'
                subtext = varargin{opt+1};
                opt = opt + 2;
            case 'elementlist'
                ellist = varargin{opt+1};
            case 'normref'
                normref = varargin{opt+1};
            otherwise
                error('Unknown field type');
        end
        opt = opt + 2;
    end
    
    % set order of arcs by SortField
    [~,ia] = sort(arcs{:,sortfield}(:,1));
    ic = ia([1:height(arcs)]');
    cmap = flipud([flipud(autumn(ceil(length(ic)/2))); ...
        winter(ceil(length(ic)/2))]);
end

% default plot colors
for i = 1:height(arcs)
    switch arcs.basementType{i}
        case 'continental'
            colour(i,:) = [1 0 0];
            leg(1) = i;
        case 'oceanic'
            colour(i,:) = [0 0 1];
            leg(3) = i;
        case 'ribbon'
            colour(i,:) = [1 0.5 0];
            leg(2) = i;
    end
end
fid = fopen('spidernorm.log','w+');
fclose(fid);
figure;
aind = data.arcid > 0;

for i = 1:length(ellist)
    data{:,[lower(ellist{i}),'_ppm_adj']} = nan([height(data),1]);
end

for i = 1:height(arcs)
    fid = fopen('spidernorm.log','a+');
    fprintf(fid,'\nComputing spidernorm for %s...\n',arcs.name{i});
    fclose(fid);
    
    subset(:,i) = data.arcid == ia(i);
    subtext{i} = arcs.name{ia(i)};
    
    %if nargin == 2
    %    set(p(i),'Color',colour(i,:));
    %else
    colour(i,:) = cmap(i,:);
    %    ltxt{i} = [arcs.name{i},' (',num2str(arcs{i,sortfield}),')'];
    %end
end

t = spider2(data, ...
        'Subset',subset, ...
        'Subtext',subtext, ...
        'Color',colour, ...
        'NormRef',normref, ...
        'TraceNorm',{reffield,refval}, ...
        'Elements',ellist, 'LogFile');

hpax([-1 3],'y');
ylabel(['Abundance/NMORB (',normref,')']);

return