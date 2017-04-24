function [values, isolated] = prune_sparse(N, values, thr)
pruneP = values>=thr;
%if ~isempty(pruneP)
values = values(pruneP) ;

isolated = setdiff(1:N, nnz(values));
fprintf('%d (%.2f%%) of %d components are isolated at threshold %f\n', ...
    length(isolated), length(isolated)/N*100, N, thr);