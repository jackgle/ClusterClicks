function [specClickTf_norm,specClickTf_norm_short] = ...
    spec_norm_diff(specClickTf,stIdx,edIdx)

minSSsection = min(specClickTf,[],2);
specClickTf_minNorm = (specClickTf - ...
  minSSsection(:,ones(1,size(specClickTf,2))));
maxSSsection = max(specClickTf_minNorm,[],2);
specClickTf_norm = specClickTf_minNorm./maxSSsection(:,ones(1,size(specClickTf_minNorm,2)));
specClickTf_norm_short = specClickTf_norm(:,stIdx:edIdx);
specClickTf_norm_short = diff(specClickTf_norm_short,1,2);
