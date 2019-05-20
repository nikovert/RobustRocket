function prop = init_rocket(string)
%Reads the values from a given table and returns a structure containing the
%relevant information.
    prop = readtable(string,'Format','%s%f');
%     motorname=char(table2array(prop(end,1)));
%     motordata = rocketmotor(motorname);
    
end