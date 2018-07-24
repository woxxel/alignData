

function para = set_paras(mouse)
  
  para = struct;    %% maybe outside?
  para.f = 15;                 %% Hz frame rate of sampling
  
  %% information about data structure in behaviour files
  %%% maybe add string to search for in here
  para.startframe = 1;
  para.nframe = 8989;
  
  para.col_t = 1;
  para.col_speed = 2;
  para.col_frame = 4;
  para.col_pos = 6;
  para.col_TTL = 8;
  
  
  %% procession of behavioural data
  para.runthres=0.5;                            %% threshold to define active / inactive periods (cm/sec)
  para.lr_min = 30;                             %% number of minimal frames considered to be a "long run"-event
  para.nbin=80;                                 %% number of divisions of the linear track
  
  para.cols = 8;
  
  if any(mouse==[72,'72NG'])
    para.cols = 9;
    para.col_TTL = 3;
  end
  
  if any(mouse==[72,884])
    para.totallength = 1600;     %% length of the linear track
    para.pos_offset = 1600;   %% offset of position values
    
    para.nSes = 10;
    
  
  elseif any(mouse==[245])
    %%% mouse 245   %% should the end and start position be contained in the "active" data? -> bias towards information during start/end position!?
    para.totallength = 600;     %% length of the linear track
    para.pos_offset = 300;   %% offset of position values
    
    para.nSes = 30;
    
  end
  
  para.binwidth = para.totallength/para.nbin;
  
end