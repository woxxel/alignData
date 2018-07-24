

function [bh] = prep_behaviour(para,pathMouse,loadData)
  
  plt = true;
  disp('preparing behavioral data...')
  
  pathSave = pathcat(pathMouse,'behaviour_data.mat');
  if exist(pathSave,'file') && loadData
    disp(sprintf('Behavior file already exists. Loading from %s...',pathSave))
    load(pathSave)
    return
  end
  
  bh = struct('duration',cell(para.nSes,1),...
              'time',zeros(1,para.nframe),...
              'frame',zeros(1,para.nframe),...
              'position',zeros(1,para.nframe),...
              'speed',zeros(1,para.nframe),...
              'runrest',zeros(1,para.nframe),...
              'longrunperiod',zeros(1,para.nframe),...
              'dwelltime',zeros(1,para.nbin),...
              'norm_dwelltime',zeros(1,para.nbin));
  
  for s=1:para.nSes
    disp(sprintf('Processing session %d',s))
    
    pathSaveSes = pathcat(pathMouse,sprintf('Session%02d',s),'behaviour_data.mat');
    
    % load data
    pathSession = pathcat(pathMouse,sprintf('Session%02d',s));
    bhfile=dir(pathcat(pathSession,'*.txt'));
    pathBH = pathcat(pathSession,bhfile.name);
    bh(s) = read_BH_file(bh(s),pathBH,para);              % reading and processing behaviour data
    bh(s).binpos = min(floor(bh(s).position/para.binwidth)+1,para.nbin);
    
    %% defining run epochs by finding super-threshold movement speed
    sm = ones(1,10+1);
    bh(s).runrest = imerode(imdilate(bh(s).speed >= para.runthres,sm),sm);   %% filling holes of size up to 5 (=1/3 sec)
    bh(s).longrunperiod = bwareaopen(bh(s).runrest,para.lr_min);    % find lonrunperiods as connected areas > lr_min of size
    
    %% create dwell time histogram
    bh(s).dwelltime = zeros(1,para.nbin);
    binpos = bh(s).binpos(bh(s).longrunperiod);
    for i=1:length(binpos)
      bh(s).dwelltime(binpos(i)) = bh(s).dwelltime(binpos(i)) + 1/para.f;
    end
    bh(s).norm_dwelltime = bh(s).dwelltime/sum(bh(s).dwelltime);
    
    
    if strcmp(para.mode,'regress')
      binpos = min(floor(bh(s).position/para.binwidth)+1,para.nbin);
      bh(s).place_vectors = zeros(para.nbin,para.nframe);
      for i = 1:length(binpos)
        if bh(s).longrunperiod(i)
          bh(s).place_vectors(binpos(i),i) = 1;
        end
      end
    end
    
    
    
%      if plt
%        close all
%        figure('position',[100 100 1500 700])
%        ymax = ceil(max(bh(s).speed));
%        
%        yyaxis left
%        hold on
%        bar(bh(s).time,bh(s).longrunperiod*ymax,1,'FaceColor',[0.8 0.8 0.8])
%        plot(bh(s).time,bh(s).position*ymax/para.totallength,'b')
%        plot([0,bh(s).duration],[60,60]*ymax/para.nbin,'r--','LineWidth',2)
%        plot([0,bh(s).duration],[20,20]*ymax/para.nbin,'g--','LineWidth',2)
%        hold off
%        ylim([-ymax ymax])
%        yticks([0,ymax])
%        yticklabels([0,1])
%        ylabel('Position')
%        
%        yyaxis right
%        plot(bh(s).time,-bh(s).speed,'r')
%        ylim([-ymax ymax])
%        yticks([-ymax,0])
%        yticklabels([-ymax,0])
%        ylabel('velocity [cm/s]')
%        xlabel('time [s]')
%        
%  %        pathName = pathcat(pathMouse,sprintf('behavior_%02d.png',s));
%  %        print(pathName,'-dpng','-r600')
%        waitforbuttonpress;
%      end
    
    crop_behavior = bh(s);
    save(pathSaveSes,'crop_behavior','-v7.3')
    disp(sprintf('Behavioural data saved under %s',pathSaveSes))
  end
  
  save(pathSave,'bh','-v7.3')
  disp(sprintf('Behavioural data saved under %s',pathSave))

end



function [bh] = read_BH_file(bh,pathBH,para)
  % This program extracts the part of image acquisition from the whole behavioral data. 
  % The existance of frame numbers is examined and if not present, new
  % numbers are calculated from the TTL signal and assigned.
  %
  % INPUT ARGUMENT
  % pathBH    - path to the behavior data file
  % para      - structure containing all kind of parameters for the analysis
  
  % OUTPUT ARGUMENT
  % bh        - structure containing data at each time frame: time, frame number, position, speed and overall duration of session
  %
  % 2018-06-07 by Alexander Schmidt
  
  %%% opening and reading data from pathBH
  
%    delimiter = sprintf('\t',''); %or whatever
%    fid = fopen(pathBH,'r');
%    tLines = fgets(fid);
%    numCols = numel(strfind(tLines,delimiter)) + 1;
%    fclose(fid);
%    numCols
  
  fid = fopen(pathBH,'r');
  fgets(fid);
  whole_data=fscanf(fid, '%f', [para.cols, inf]);
  fclose(fid);
  
  
  %%% check, whether frame numbers are present and whether they are in sync
  data_frame  = whole_data(para.col_frame,:);     % frame number
  sync_sum    = sum(data_frame);
  sync_last   = data_frame(end);
  
  data_TTL    = whole_data(para.col_TTL,:);       % microscope TTL
  
  if sync_sum~=0 & para.nframe==sync_last
      
      idx_start=find(data_frame==para.startframe,1,'first');
      idx_end=find(data_frame==para.nframe-1,1,'last')+3;
      
  elseif  sync_sum==0 | para.nframe~=sync_last
      
      idx_start=find(data_TTL<0.5,1,'first');% find the first point or low TTL    
      diffsig=diff(data_TTL);
      idx_end=find(diffsig > 0.5,1,'last');% find the last point or high TTL
            
      if sync_sum==0
              disp('WARNING --- Sync signal is missing')
      elseif para.nframe~=sync_last
              disp('WARNING --- Frame number mismatch')
      end
  end
  
  datapoints = idx_end-idx_start+1;
  %% if needed, reassign frame numbers
  if  sync_sum==0 | para.nframe~=sync_last
    datapoints_per_frame  = datapoints/para.nframe; 
    frame_num = ceil((1:datapoints)/datapoints_per_frame);
    disp('New frame numbers assigned')
  else
    frame_num = whole_data(para.col_frame,idx_start:idx_end);
  end
  
%    whole_data(para.col_t,1:100)
  %%% read remaining data from file
  bh_raw = struct;
  bh_raw.datapoints = datapoints;
  bh_raw.time = whole_data(para.col_t,idx_start:end)-whole_data(para.col_t,idx_start);       % offset corrected
  bh_raw.pos = whole_data(para.col_pos,idx_start:end) + para.pos_offset;
  bh_raw.speed = whole_data(para.col_speed,idx_start:end);
  
  align = true;
  if align
    setappdata(0,'behavior_data',bh_raw)
    setappdata(0,'parameter',para)
    
    uiwait(alignData)
    disp('done')
    bh_raw.pos = getappdata(0,'new_position');
    rmappdata(0,'behavior_data')
    rmappdata(0,'parameter')
    rmappdata(0,'new_position')
  end
  
  bh_raw.time = bh_raw.time(1:datapoints);
  bh_raw.pos = bh_raw.pos(1:datapoints);
  bh_raw.speed = bh_raw.speed(1:datapoints);
  
  bh_raw.speed = imgaussfilt(bh_raw.speed,10);                                                    % filter speed with 200ms window
  
  
  %% resample stuff at 15 Hz by averaging
  for i=1:para.nframe
    idx = find(frame_num==i);
    bh.time(i) = mean(bh_raw.time(idx));        % time
    bh.frame(i) = i;                            % frame number
    bh.position(i) = median(bh_raw.pos(idx));   % position
    bh.speed(i) = mean(bh_raw.speed(idx));      % speed
  end
  bh.duration = bh.time(end);
  
end