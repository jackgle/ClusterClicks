function [values, rows, cols, isolated] = prune(N, values, rows, cols, thr)
pruneP = values>=thr;
%if ~isempty(pruneP)
values = values(pruneP) ;
rows =  rows(pruneP);
cols =  cols(pruneP);
%end

isolated = setdiff(1:N, unique([rows; cols]));
fprintf('%d (%.2f%%) of %d components are isolated at threshold %f\n', ...
    length(isolated), length(isolated)/N*100, N, thr);