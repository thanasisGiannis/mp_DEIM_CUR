function [A] = load_mat_ex_2()
filename = 'GSE10072_series_matrix.txt';
[A,~,~] = load_GSE10072(filename);
A = A - mean(A,2);
end

function [A, probeIDs, sampleIDs] = load_GSE10072(filename)
%LOAD_GSE10072 Load the GSE10072 GEO series matrix.
%
%   [A, probeIDs, sampleIDs] = load_GSE10072(filename)
%
%   A          : 22283 x 107 expression matrix, row-centered
%   probeIDs   : probe identifiers
%   sampleIDs  : GEO sample identifiers

lines = readlines(filename);

% Locate the expression table
ibegin = find(strtrim(lines) == "!series_matrix_table_begin", 1);
iend   = find(strtrim(lines) == "!series_matrix_table_end", 1);

if isempty(ibegin) || isempty(iend)
    error('Could not find GEO series matrix table.');
end

% Header contains ID_REF followed by sample identifiers
header = split(lines(ibegin + 1), sprintf('\t'));
header = erase(header, '"');

sampleIDs = header(2:end);

% Data lines
dataLines = lines(ibegin + 2 : iend - 1);

m = length(dataLines);
n = length(sampleIDs);

A = zeros(m,n);
probeIDs = strings(m,1);

for i = 1:m
    fields = split(dataLines(i), sprintf('\t'));
    fields = erase(fields, '"');

    probeIDs(i) = fields(1);
    A(i,:) = str2double(fields(2:end)).';
end

% Sorensen--Embree Example 3:
% subtract the mean of each row
A = A - mean(A,2);

fprintf('Loaded matrix: %d x %d\n', size(A,1), size(A,2));

end

