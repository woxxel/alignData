
function check_bh(basePath,mouse,thr,nsd,mode)

  loadData = false;
  
  pathData = sprintf('%s%d',basePath,mouse)
%    loadPath = pathcat(pathData,sprintf('ROI_final_p=%4.2g.mat',p_thr))
%    load(loadPath)
%    suffix = sprintf('_s=%d_p=%4.2g',thr(1),thr(2));
  
  
  para = set_paras(mouse,mode,thr,nsd);
  
  bh = prep_behaviour(para,pathData,loadData);

end