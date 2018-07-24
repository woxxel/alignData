
function check_bh(basePath,mouse)

  loadData = false;
  
  pathData = sprintf('%s%s',basePath,string(mouse))
  para = set_paras(mouse);
  bh = prep_behaviour(para,pathData,loadData);

end