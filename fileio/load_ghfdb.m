function data = load_ghfdb

%fmt = '%f%f%f%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%f%s%f%f%f%f%f%s%s%f%f%f%s%s%s%f%s%f%f%s%s%s%s%s%s%s%s%s%s';
data = readtable('../database/heat_flow/master.xlsx');%,'Format',fmt);

return