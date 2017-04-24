function [values, rows, cols, isolated] = prune_by_node1(N, values, rows, cols, thr)
pRe = squareform(values);
nPerc = .05;
nLinks = 100;%round(nPerc*N);
fprintf('Retaining top %d%% linkages per node = %d linkages\n',round(nPerc*100),nLinks)

for iR = 1:length(pRe)
    [~,I] = sort(pRe(iR,:),'descend');
    pRe(iR,I(nLinks+1:end)) = 0;
    % fprintf('done with row %d \n',iR)
end
values = squareform(pRe);
pruneP = values>=thr;
%if ~isempty(pruneP)
values = values(pruneP) ;
rows =  rows(pruneP);
cols =  cols(pruneP);
% end

isolated = setdiff(1:N, unique([rows; cols]));
fprintf('%d (%.2f%%) of %d components are isolated at threshold %f\n', ...
    length(isolated), length(isolated)/N*100, N, thr);