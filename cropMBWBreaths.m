function [ table ] = cropMBWBreaths(table,  mbwBreaths, numberOfBreaths )
%Takes a structure, iterates through its fields and cuts the breaths to the
%right length
%   Detailed explanation goes here
    fields=fieldnames(table);
    for i=1:numel(fields)
        if length(table.(fields{i}))==numberOfBreaths
            table.(fields{i})=table.(fields{i})(mbwBreaths,:);
        end
        if isa(table.(fields{i}),'Gas')
            table.(fields{i})=cropMBWBreaths(table.(fields{i}), mbwBreaths, numberOfBreaths);
        end
    end
end

