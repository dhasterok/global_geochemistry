function varargout = ts2age(varargin)

if nargin > 1
    error('There may be at most one arguement to ts2age.');
end

%fid = fopen('geologic_time_scale.csv');
%ts = textscan(fid,'%s %d %d','Delimiter',',');
%fclose(fid);

ts = timescale;

if nargin == 0
    varargout = {ts};
    return
else
    for i = 1:length(ts(:,1))
        ts{i,1} = lower(ts{i,1});
    end
    agename = varargin{1};
    if length(agename) < 1
        ages = [NaN NaN NaN];
    end
end

if ischar(agename)
    ages = findage(agename,ts);
elseif iscell(agename)
    for i = 1:length(agename)
        ages = findage(agename{i},ts);
    end
end

varargout = {ages};

return

function age = findage(agename,ts)

i = strmatch(lower(agename),ts(:,1));
if isempty(i)
    age = [NaN NaN NaN];
else
    age = [ts{i,2} 0.5*(ts{i,2}+ts{i,3}) ts{i,3}];
end

return



function ts = timescale;


ts = {'Phanerozoic',0,541;
    'Cenozoic',0,66;
    'Quaternary',0,2.6;
    'Holocene',0,0.01;
    'Pleistocene',0.01,2.6;
    'Calabrian',0.01,1.8;
    'Gelasian',1.8,2.6;
    'Neogene',2.6,23;
    'Pliocence',2.6,5.3;
    'Piacenzian',2.6,3.6;
    'Zanclean',3.6,5.3;
    'Miocene',5.3,23;
    'Messinian',5.3,7.2;
    'Tortonian',7.2,11.6;
    'Serravallian',11.6,13.8;
    'Langhian',13.8,16;
    'Burdigalian',16,20.4;
    'Aquitanian',20.4,23;
    'Paleogene',23,66;
    'Oligocene',23,33.9;
    'Chattian',23,28.1;
    'Rupelian',28.1,33.9;
    'Eocene',33.9,56;
    'Priabonian',33.9,37.8;
    'Bartonian',37.8,41.2;
    'Lutetian',41.2,47.8;
    'Ypresian',47.8,56;
    'Paleocene',56,66;
    'Thanetian',56,59.2;
    'Selandian',59.2,61.6;
    'Danian',61.6,66;
    'Mesozoic',66,252;
    'Late Cretaceous',66,100;
    'Maastrichtian',66,72.1;
    'Campanian',72.1,83.6;
    'Santonian',83.6,86.3;
    'Coniacian',86.3,89.8;
    'Turonian',89.8,93.9;
    'Cenomanian',93.9,100;
    'Early Cretaceous',100,145;
    'Albian',100,113;
    'Aptian',113,126;
    'Barremian',126,131;
    'Jurassic',145,201;
    'Late Jurassic',145,164;
    'Hauterivian',131,134;
    'Valanginian',134,139;
    'Berriasian',139,145;
    'Tithonian',145,152;
    'Kimmeridgian',152,157;
    'Oxfordian',157,164;
    'Middle Jurassic',164,174;
    'Callovian',164,166;
    'Bathonian',166,168;
    'Bajocian',168,170;
    'Aalenian',170,174;
    'Early Jurassic',174,201;
    'Toarcian',174,183;
    'Pliensbachian',183,191;
    'Sinemurian',191,199;
    'Hettangian',199,201;
    'Triassic',201,252;
    'Late Triassic',201,237;
    'Rhaetian',201,209;
    'Norian',209,228;
    'Carnian',228,237;
    'Middle Triassic',237,247;
    'Ladinian',237,241;
    'Anisian',241,247;
    'Early Triassic',247,252;
    'Olenekian',247,250;
    'Induan',250,252;
    'Paleozoic',252,541;
    'Permian',252,299;
    'Lopingian',252,260;
    'Changhsingian',252,254;
    'Wuchiapingian',254,260;
    'Guadalupian',260,272;
    'Capitanian',260,265;
    'Wordian',265,269;
    'Roadian',269,272;
    'Cisuralian',272,299;
    'Kungurian',272,279;
    'Artinskian',279,290;
    'Sakmarian',290,296;
    'Asselian',296,299;
    'Carboniferous',299,359;
    'Pennsylvanian',299,323;
    'Late Pennsylvanian',299,307;
    'Gzhelian',299,304;
    'Kasimovian',307,307;
    'Middle Pennsylvanian',307,315;
    'Moscovian',307,315;
    'Early Pennsylvanian',315,323;
    'Bashkirian',315,323;
    'Mississippian',323,359;
    'Late Mississippian',323,331;
    'Serpukhovian',323,331;
    'Middle Mississippian',331,347;
    'Visean',331,347;
    'Early Mississippian',347,359;
    'Tournaisian',347,359;
    'Devonian',359,419;
    'Late Devonian',359,383;
    'Famennian',359,372;
    'Frasnian',372,383;
    'Middle Devonian',383,393;
    'Givetian',383,388;
    'Eifelian',388,393;
    'Early Devonian',393,419;
    'Emsian',393,408;
    'Pragian',408,411;
    'Lochkovian',411,419;
    'Silurian',419,444;
    'Pridoli',419,423;
    'Ludlow',423,427;
    'Ludfordian',423,426;
    'Gorstian',426,427;
    'Wenlock',430,430;
    'Homerian',427,430;
    'Sheinwoodian',430,433;
    'Llandovery',433,444;
    'Telychian',433,439;
    'Aeronian',439,441;
    'Rhuddanian',441,444;
    'Ordovician',444,485;
    'Late Ordivician',444,458;
    'Hirnantian',444,445;
    'Katian',445,453;
    'Sandbian',453,458;
    'Middle Ordivician',458,470;
    'Darriwilian',458,467;
    'Dapingian',467,470;
    'Early Ordivician',470,485;
    'Floian',470,478;
    'Tremadocian',478,485;
    'Cambrian',485,541;
    'Furongian',485,497;
    'Age 10',485,490;
    'Jiangshanian',490,494;
    'Paibian',494,497;
    'Epoch 3',497,509;
    'Guzhangian',497,501;
    'Drumian',501,505;
    'Age 5',505,509;
    'Epoch 2',509,521;
    'Age 4',509,514;
    'Age 3',514,521;
    'Terreneuvian',521,541;
    'Age 2',521,529;
    'Fortunian',529,541;
    'Precambrian',541,4000;
    'Proterozoic',541,2500;
    'Neoproterozoic',541,1000;
    'Ediacaran',541,635;
    'Cryogenian',635,860;
    'Tonian',860,1000;
    'Mesoproterozoic',1000,1600;
    'Stenian',1000,1200;
    'Ectasian',1200,1400;
    'Calymmian',1400,1600;
    'Paleoproterozoic',1600,2500;
    'Statherian',1600,1800;
    'Orosirian',1800,2050;
    'Rhyacian',2050,2300;
    'Siderian',2300,2500;
    'Archean',2500,4000;
    'Neoarchean',2500,2800;
    'Mesoarchean',2800,3600;
    'Paleoarchean',3600,4000;
    'Eoarchean',4000,4650;
    'Hadean',4000,4650};

return
