%% GUI to align position and speed data
%%
%% INPUT via getappdata from "prep_behaviour.m": preloaded data from behavior-text-file
%% OUTPUT via setappdata to "prep_behavior.m": aligned position values

%% written 24.07.2018 by A. Schmidt
%%
%% TODOs:
%% - give estimate, when data seems to be alright (from overall number of alignment points?
%% - adapt y-axis to fit either bin-number of cm


function varargout = alignData(varargin)
% POS_SPEED_ALIGN MATLAB code for alignData.fig
%      POS_SPEED_ALIGN, by itself, creates a new POS_SPEED_ALIGN or raises the existing
%      singleton*.
%
%      H = POS_SPEED_ALIGN returns the handle to a new POS_SPEED_ALIGN or the handle to
%      the existing singleton*.
%
%      POS_SPEED_ALIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POS_SPEED_ALIGN.M with the given input arguments.
%
%      POS_SPEED_ALIGN('Property','Value',...) creates a new POS_SPEED_ALIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alignData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alignData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alignData

% Last Modified by GUIDE v2.5 29-Jul-2018 12:57:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @alignData_OpeningFcn, ...
                   'gui_OutputFcn',  @alignData_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before alignData is made visible.
function alignData_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alignData (see VARARGIN)
  
  h.save_status = false;
  h.data = struct;
  
  assert(nargin<=4,'Too many inputs provided. Call alignData with a mouse folder- or file-name, or without any input!')
  
  h = set_paras(h);
  
  if nargin == 4   %% if GUI is called with input
    
    %% check if input is string
    h.pathIn = varargin{1};
    assert(isstr(h.pathIn),'Provide input as a string of the data path')
    
    h.currSession = 1;
    
    [folderName,fileName,extension] = fileparts(h.pathIn);
    if strcmp(extension,'.txt')
      h.pathSession = folderName;
      h = read_bh_file(h,h.pathIn);
      h.nSessions = 0;
    else
      h.folders = dir(pathcat(h.pathIn,'Session*'));
      h.nSessions = length(h.folders);
      h.pathSession = pathcat(h.pathIn,h.folders(h.currSession).name);
      
      bhfile=dir(pathcat(h.pathIn,h.folders(h.currSession).name,'*m.txt'));
      h.pathBHfile = pathcat(h.pathIn,h.folders(h.currSession).name,bhfile.name);
      h = read_bh_file(h,h.pathBHfile);
      set(h.button_next,'string','>>')
      set(h.button_previous,'visible','on')
    end
    h = run_methods(h);
    
  else
    h.nSessions = 0;
    h.currSession = 0;
    h.pathSession = '';
  end
  
%    set(h.entry_col_reward,'enable','off')
  
  % Choose default command line output for alignData
  h.output = hObject;

  % Update h structure
  guidata(hObject, h);
  
% UIWAIT makes alignData wait for user response (see UIRESUME)
% uiwait(h.alignData);


function show_info()
  
  msgbox(sprintf('Info for usage\n\nStop onsets (anchors) in speed (red dots) and position (black dots) arrays are obtained by finding periods of 0 value (speed) or 0 change (position) in smoothed data with a gaussian kernel of length "SMOOTH"\n\nClick on anchor points to remove / add them to the alignment procedure. The updated best guess of aligned data (black line) will be displayed immediately'))
  

function h = set_paras(h)
  
  h.paras = struct;    %% maybe outside?
  h.paras.f = 15;                 %% Hz frame rate of sampling
  
  %% information about data structure in behaviour files
  %%% maybe add string to search for in here
  h.paras.startframe = 1;
  h.paras.nframe = 8989;
  
  %% structure of txt-file -> dynamic?
  h.paras.col_t = 1;
  h.paras.col_speed = 2;
  h.paras.col_frame = 4;
  h.paras.col_pos = 6;
  
  if isempty(get(h.entry_col_reward,'string'))
    h.paras.col_reward = 8;   %% -> dynamic
    set(h.entry_col_reward,'string',sprintf('%d',h.paras.col_reward))
  else
    h.paras.col_reward = str2double(get(h.entry_col_reward,'string'));
    
  end
  
  if isempty(get(h.entry_col_TTL,'string'))
    h.paras.col_TTL = 3;
    set(h.entry_col_TTL,'string',sprintf('%d',h.paras.col_TTL))
  else
    h.paras.col_TTL = str2double(get(h.entry_col_TTL,'string'));
  end
  
  
  %% procession of behavioural data
  h.paras.runthres = 0.5;                            %% threshold to define active / inactive periods (cm/sec)
  h.paras.lr_min = 30;                             %% number of minimal frames considered to be a "long run"-event
  h.paras.nbin = 80;                                 %% number of divisions of the linear track
  
  h.paras.cols = 8;
  
  h.paras.cmlength = 100;
  
  
  
function h = set_paras_from_entries(h)
  h.paras.ptlength = str2double(get(h.entry_tracklength,'string'));
  h.paras.pos_offset = str2double(get(h.entry_pos_offset,'string'));
  
  if isempty(get(h.entry_reward_position,'string'))
    h.paras.reward_thr = 3/4*h.paras.ptlength;
    set(h.entry_reward_position,'string',h.paras.reward_thr/h.paras.ptlength*h.paras.cmlength)
  else
    h.paras.reward_thr = str2double(get(h.entry_reward_position,'string'))/h.paras.cmlength*h.paras.ptlength;
  end
  
  h.paras.col_reward = str2double(get(h.entry_col_reward,'string'));
  h.paras.col_TTL = str2double(get(h.entry_col_TTL,'string'));
  
  
function h = read_bh_file(h,pathBH)
  
  delimiter = sprintf('\t',''); %or whatever
  fid = fopen(pathBH,'r');
  fgets(fid);
  tLines = fgets(fid);
  numCols = numel(strfind(tLines,delimiter));
  whole_data=fscanf(fid, '%f', [numCols, inf]);
  fclose(fid);
  
  %%% check, whether frame numbers are present and whether they are in sync
  data_frame  = whole_data(h.paras.col_frame,:);     % frame number
  sync_sum    = sum(data_frame);
  sync_last   = data_frame(end);
  data_reward = whole_data(h.paras.col_TTL,:);       % microscope TTL
  
  if sync_sum~=0 & h.paras.nframe==sync_last
      
      idx_start=find(data_frame==h.paras.startframe,1,'first');
      idx_end=find(data_frame==h.paras.nframe-1,1,'last')+3;
      
  elseif  sync_sum==0 | h.paras.nframe~=sync_last
      
      idx_start=find(data_reward<0.5,1,'first');% find the first point or low reward   
      diffsig=diff(data_reward);
      idx_end=find(diffsig > 0.5,1,'last');% find the last point or high reward
  end
  
  datapoints = idx_end-idx_start+1;
  %% if needed, reassign frame numbers
  if  sync_sum==0 | h.paras.nframe~=sync_last
    datapoints_per_frame  = datapoints/h.paras.nframe; 
    h.data.frame_num = ceil((1:datapoints)/datapoints_per_frame);
  else
    h.data.frame_num = whole_data(h.paras.col_frame,idx_start:idx_end);
  end
  
  %% raw data
  h.data.bh.datapoints = datapoints;
  h.data.bh.time = whole_data(h.paras.col_t,idx_start:end)-whole_data(h.paras.col_t,idx_start);       % offset corrected
  h.data.bh.speed = whole_data(h.paras.col_speed,idx_start:end);
  h.data.bh.pos = whole_data(h.paras.col_pos,idx_start:end);
  h.data.bh.TTL = whole_data(h.paras.col_TTL,idx_start:end);
  h.data.bh.reward = whole_data(h.paras.col_reward,idx_start:end);
  
  min_pos = min(h.data.bh.pos);
  max_pos = max(h.data.bh.pos);
  
  if isempty(get(h.entry_tracklength,'string'))
    h.paras.ptlength = round((max_pos-min_pos)/100)*100;
    set(h.entry_tracklength,'string',sprintf('%d',h.paras.ptlength))
  end
  if isempty(get(h.entry_pos_offset,'string'))
    h.paras.pos_offset = -round(min_pos);
    set(h.entry_pos_offset,'string',sprintf('%d',h.paras.pos_offset))
  end
  
  %% cropped, aligned data
  h.data.time = h.data.bh.time(1:h.data.bh.datapoints);
  h.data.duration = h.data.time(end);
  h.data.speed = h.data.bh.speed(1:h.data.bh.datapoints);  
  
  
  
function h = run_methods(h)
  
  h = set_paras_from_entries(h);
  h = smooth_data(h);
  h = run_align(h);
  h = plot_align(h);
  h = updateSaveStatus(h,false);
  

function h = smooth_data(h)
  
  % smooth data to obtain longer stopping periods, only
  smooth_para = str2double(get(h.entry_smooth,'string'))*(h.data.bh.datapoints/h.data.duration);
  pos_smooth = imgaussfilt(h.data.bh.pos,smooth_para);
  speed_smooth = imgaussfilt(h.data.bh.speed,smooth_para);
  
  %%% find resting periods in both, speed and location vector
  dx = diff(pos_smooth);
  dx_bool = (dx == 0);
  idx_run_pos = find(~dx_bool);
  h.data.idx_align_pos = [0,idx_run_pos(diff(idx_run_pos)>1)];
  h.data.idx_align_pos_mask = true(size(h.data.idx_align_pos));
  
  h.data.speed_bool = (speed_smooth == 0);
  idx_run_speed = find(~h.data.speed_bool)+1;
  h.data.idx_align_speed = [0,idx_run_speed(diff(idx_run_speed)>1)];
  h.data.idx_align_speed_mask = true(size(h.data.idx_align_speed));
  
  reward_signal = find(bwareaopen(h.data.bh.reward>0.5,2));
  h.data.idx_reward_signal = [0, reward_signal(find(diff(reward_signal)>1))];
  
  reward_pos = find((h.data.bh.pos+h.paras.pos_offset)<h.paras.reward_thr);
  h.data.idx_reward_pos = [0, reward_pos(find(diff(reward_pos)>1))];
  
  
function h = run_align(h)
  
  h.data.pos_s = zeros(1,length(h.data.bh.pos));
  offset = 1;
  i = 2;
  
  idx_align_pos = h.data.idx_align_pos(h.data.idx_align_pos_mask);
  idx_align_speed = h.data.idx_align_speed(h.data.idx_align_speed_mask);
  
  while offset < h.data.bh.datapoints
    n_loc = idx_align_pos(i)-idx_align_pos(i-1);
    n_speed = idx_align_speed(i)-idx_align_speed(i-1);
    %% check if distances of stopping periods are the same
    
    if n_loc~=n_speed   %% if not, resample location points and write to new position
      loc_tmp = h.data.bh.pos(idx_align_pos(i-1)+1:idx_align_pos(i)) + h.paras.pos_offset;
      
      trialpoint = [0,find(diff(loc_tmp)<-10),n_loc];
      for t = 2:length(trialpoint)
        
        loc_tmp_tp = loc_tmp(trialpoint(t-1)+1:trialpoint(t));
        
        
        dloc = diff(loc_tmp_tp);
        if ~all(dloc>=0)
          disp('interp pre')
          loc_tmp
          dloc
        end
        
        idx1 = round(idx_align_speed(i-1)+1 + trialpoint(t-1)*(n_speed/n_loc));
        idx2 = round(idx_align_speed(i-1)+1 + (trialpoint(t)-1)*(n_speed/n_loc));
        
        n_points = trialpoint(t) - trialpoint(t-1);
        n_fill = idx2-idx1+1;
        
        time_speed = linspace(h.data.bh.time(idx1),h.data.bh.time(idx2),n_fill);
        time_loc = linspace(h.data.bh.time(idx1),h.data.bh.time(idx2),n_points);
        
        h.data.pos_s(offset:offset+n_fill-1) = interp1(time_loc,loc_tmp_tp,time_speed);
        offset = offset + n_fill;
        
        fill_tmp = interp1(time_loc,loc_tmp_tp,time_speed);
        dloc = diff(fill_tmp);
        if ~all(dloc>=0)
          disp('interp post')
          loc_tmp
          dloc
        end
        
      end
    
    else      %% if they are, just write old to new data
      loc_tmp = h.data.bh.pos(idx_align_pos(i-1)+1:idx_align_pos(i)) + h.paras.pos_offset;
      
      trialpoint = [0,find(diff(loc_tmp)<-10),n_loc];
      for t = 2:length(trialpoint)
        
        loc_tmp_tp = loc_tmp(trialpoint(t-1)+1:trialpoint(t));
        dloc = diff(loc_tmp_tp);
        if ~all(dloc>=0)
          disp('direct')
          loc_tmp
          dloc
        end
        n_points = trialpoint(t) - trialpoint(t-1);
        h.data.pos_s(offset:offset+n_points-1) = loc_tmp_tp;
        offset = offset + n_points;
        
      end
      
    end
    
    %% if there is no stopping period after recording stops, assign remaining positions linearly
    if i == length(idx_align_pos) || i == length(idx_align_speed) || offset > h.data.bh.datapoints
      lbh = length(h.data.bh.pos);
      n_speed = min(lbh - idx_align_speed(i),lbh - idx_align_pos(i));
      
      if n_speed > 0
        h.data.pos_s(offset:offset+n_speed-1) = h.data.bh.pos(idx_align_pos(i)+1:idx_align_pos(i)+n_speed) + h.paras.pos_offset;
      end
      break
    end
    i = i+1;
    
  end
  
  h.data.pos_r = zeros(1,length(h.data.bh.pos));
  
  if length(h.data.idx_reward_signal)>1 && length(h.data.idx_reward_pos)>1
    
    offset = 1;
    i = 2;
    while offset < h.data.bh.datapoints
      
      n_pos = h.data.idx_reward_pos(i)-h.data.idx_reward_pos(i-1);
      n_signal = h.data.idx_reward_signal(i)-h.data.idx_reward_signal(i-1);
      
      if n_pos~=n_signal   %% if not, resample location points and write to new position
        
        loc_tmp = h.data.bh.pos(h.data.idx_reward_pos(i-1)+1:h.data.idx_reward_pos(i)) + h.paras.pos_offset;
        
        trialpoint = [0,find(diff(loc_tmp)<-10),n_pos];
        for t = 2:length(trialpoint)
          loc_tmp_tp = loc_tmp(trialpoint(t-1)+1:trialpoint(t));
          
          idx1 = round(h.data.idx_reward_signal(i-1)+1 + trialpoint(t-1)*(n_signal/n_pos));
          idx2 = round(h.data.idx_reward_signal(i-1)+1 + (trialpoint(t)-1)*(n_signal/n_pos));
          
          n_points = trialpoint(t) - trialpoint(t-1);
          n_fill = idx2-idx1+1;
          
          time_signal = linspace(h.data.bh.time(idx1),h.data.bh.time(idx2),n_fill);
          time_loc = linspace(h.data.bh.time(idx1),h.data.bh.time(idx2),n_points);
          
          h.data.pos_r(offset:offset+n_fill-1) = interp1(time_loc,loc_tmp_tp,time_signal);
          offset = offset + n_fill;
          
        end
      
      else      %% if they are, just write old to new data
        loc_tmp = h.data.bh.pos(h.data.idx_reward_pos(i-1)+1:h.data.idx_reward_pos(i)) + h.paras.pos_offset;
        
        trialpoint = [0,find(diff(loc_tmp)<-10),n_pos];
        for t = 2:length(trialpoint)
          loc_tmp_tp = loc_tmp(trialpoint(t-1)+1:trialpoint(t));
          
          n_points = trialpoint(t) - trialpoint(t-1);
          h.data.pos_r(offset:offset+n_points-1) = loc_tmp_tp;
          offset = offset + n_points;
          
        end
        
%          h.data.pos_r(offset:offset+n_signal-1) = loc_tmp;
        
      end
      
      if i == length(h.data.idx_reward_pos) || i == length(h.data.idx_reward_signal) || offset > h.data.bh.datapoints
        lbh = length(h.data.bh.pos);
        n_signal = min(lbh - h.data.idx_reward_signal(i),lbh - h.data.idx_reward_pos(i));
        
        if n_signal > 0
          h.data.pos_r(offset:offset+n_signal-1) = h.data.bh.pos(h.data.idx_reward_pos(i)+1:h.data.idx_reward_pos(i)+n_signal) + h.paras.pos_offset;
        end
        break
      end
      i = i+1;
    end
  end
  
  
function h = plot_align(h)
  
  cla(h.ax,'reset')
  cla(h.ax2,'reset')
  cla(h.ax3,'reset')
  
  h.data.ymax = ceil(max(h.data.speed));
  
  plot(h.ax,h.data.time,-imgaussfilt(h.data.speed,10),'b','Hittest','off')
  
  
  
  hold(h.ax,'on')
  bar(h.ax,h.data.bh.time,~h.data.speed_bool*3*h.data.ymax,1,'FaceColor',[0.8 0.8 0.8],'Hittest','off')
  
  if get(h.radio_stops,'Value')
    LW_s = 2;
    LS_s = '-';
    anchor_vis = 'on';
    err_s_vis = 'on';
    LW_r = 0.5;
    LS_r = '--';
    err_r_vis = 'off';
  else
    LW_s = 0.5;
    LS_s = '--';
    err_s_vis = 'off';
    anchor_vis = 'off';
    LW_r = 2;
    LS_r = '-';
    err_r_vis = 'on';
  end
  
  pos_tmp = h.data.pos_s*3*h.data.ymax/h.paras.ptlength;
  h.pos = plot(h.ax,h.data.bh.time,pos_tmp,'k','Hittest','off','DisplayName','animal position (aligned)','LineStyle',LS_s,'LineWidth',LW_s);
  
  pos_r_tmp = h.data.pos_r*3*h.data.ymax/h.paras.ptlength;
  h.pos_r = plot(h.ax,h.data.bh.time,pos_r_tmp,'k','Hittest','off','DisplayName','animal position (aligned, TTL)','LineStyle',LS_r,'LineWidth',LW_r);
  
%    h.pos0 = plot(h.ax,h.data.time,(h.data.bh.pos(1:h.data.bh.datapoints) + h.paras.pos_offset)*3*h.data.ymax/h.paras.ptlength,'r--','Hittest','off','DisplayName','animal position');
  h.pos0 = plot(h.ax,h.data.bh.time,(h.data.bh.pos + h.paras.pos_offset)*3*h.data.ymax/h.paras.ptlength,'r--','Hittest','off','DisplayName','animal position');
  
  
  h.reward = plot(h.ax,[0,h.data.duration],[h.paras.reward_thr,h.paras.reward_thr]*3*h.data.ymax/h.paras.ptlength,'g--','LineWidth',2,'Hittest','off','DisplayName','reward position');
  h.gate = plot(h.ax,[0,h.data.duration],[20,20]*3*h.data.ymax/h.paras.nbin,'r--','LineWidth',2,'Hittest','off','DisplayName','gate position');
  
  plot(h.ax,[h.data.duration, h.data.duration],[-h.data.ymax,3*h.data.ymax],'k--')
  
  h.scatter_pos0 = scatter(h.ax,h.data.bh.time(h.data.idx_align_pos+1),ones(1,length(h.data.idx_align_pos)),'ko','Hittest','off','visible',anchor_vis);
  h.scatter_speed0 = scatter(h.ax,h.data.bh.time(h.data.idx_align_speed+1),-ones(1,length(h.data.idx_align_speed)),'ro','Hittest','off','visible',anchor_vis);
  
  h.scatter_pos = scatter(h.ax,h.data.bh.time(h.data.idx_align_pos(h.data.idx_align_pos_mask)+1),ones(1,sum(h.data.idx_align_pos_mask)),'ko','filled','Hittest','off','DisplayName','stops (position)','visible',anchor_vis);
  h.scatter_speed = scatter(h.ax,h.data.bh.time(h.data.idx_align_speed(h.data.idx_align_speed_mask)+1),-ones(1,sum(h.data.idx_align_speed_mask)),'ro','filled','Hittest','off','DisplayName','stops (speed)','visible',anchor_vis);
  
  hold(h.ax,'off')
  ylim(h.ax,[-h.data.ymax 3*h.data.ymax])
  yticks(h.ax,[-h.data.ymax,0,h.data.ymax])
  yticklabels(h.ax,[h.data.ymax,0,1])
  ylabel(h.ax,'Position')
  xlim(h.ax,[0,min(600,h.data.bh.time(end))])
  
  set(h.ax,'ButtonDownFcn',{@remove_align_point,h.ax},'Hittest','on','PickableParts','All');
  legend(h.ax,[h.pos0,h.pos,h.reward,h.gate,h.scatter_pos,h.scatter_speed],'location','SouthWest')
  
  cum_pos = zeros(size(h.data.bh.time));
  cum_pos(h.data.idx_align_pos(h.data.idx_align_pos_mask)+1) = 1;
  cum_pos = cumsum(cum_pos);
  
  cum_speed = zeros(size(h.data.bh.time));
  cum_speed(h.data.idx_align_speed(h.data.idx_align_speed_mask)+1) = 1;
  cum_speed = cumsum(cum_speed);
  
  cum_diff = cum_pos-cum_speed;
  
  hold(h.ax2,'on')
  plot(h.ax2,[0,h.data.bh.time(end)],[0,0],'k:')
  if get(h.radio_stops,'Value')
    h.ax2_err = plot(h.ax2,h.data.bh.time,cum_diff,'r');
    plot(h.ax2,[h.data.duration, h.data.duration],[-5,5],'k--')
    ylim(h.ax2,[-5,5])
    ylabel(h.ax2,'diff. anchor points')
  else
    err_r_0 = h.data.pos_r(1:h.data.bh.datapoints)-h.data.pos_s(1:h.data.bh.datapoints);
    err_r_0 = min([abs(err_r_0);abs(err_r_0-h.paras.pos_offset)])/(h.paras.pos_offset/100);
    h.ax2_err = plot(h.ax2,h.data.time,err_r_0,'r');
    err_ylim = max(10,max(imgaussfilt(err_r_0,10)));
    ylim(h.ax2,[-err_ylim,err_ylim])
    ylabel(h.ax2,'TTL vs stop')
  end
  
  hold(h.ax2,'off')
  xlim(h.ax2,[0,min(600,h.data.bh.time(end))])
  
  err = h.data.pos_s(1:h.data.bh.datapoints)-(h.data.bh.pos(1:h.data.bh.datapoints) + h.paras.pos_offset);
  err = min([abs(err);abs(err-h.paras.pos_offset)])/(h.paras.pos_offset/100);
  
  err_r = h.data.pos_r(1:h.data.bh.datapoints)-(h.data.bh.pos(1:h.data.bh.datapoints) + h.paras.pos_offset);
  err_r = min([abs(err_r);abs(err_r-h.paras.pos_offset)])/(h.paras.pos_offset/100);
  
  hold(h.ax3,'on')
  plot(h.ax3,[0,h.data.bh.time(end)],[0,0],'k:')
  h.err = plot(h.ax3,h.data.time,err,'k','LineStyle',LS_s,'LineWidth',LW_s);
  h.err_r = plot(h.ax3,h.data.time,err_r,'k','LineStyle',LS_r,'LineWidth',LW_r);
  
  err_ylim = max(10,max(imgaussfilt(err,10)));
  plot(h.ax3,[h.data.duration, h.data.duration],[-err_ylim,err_ylim],'k--')
  hold(h.ax3,'off')
  xlim(h.ax3,[0,min(600,h.data.bh.time(end))])
  
  ylim(h.ax3,[-err_ylim,err_ylim])
  xlabel(h.ax3,'time [s]')
  ylabel(h.ax3,'deviation [cm]')
  
  linkaxes([h.ax,h.ax2,h.ax3],'x')
  
  title(h.ax,sprintf('Now processing: %s',h.pathSession))
  
  
function update_plot(h)
  
  if get(h.radio_stops,'Value')
    LW_s = 2;
    LS_s = '-';
    anchor_vis = 'on';
    LW_r = 0.5;
    LS_r = '--';
  else
    LW_s = 0.5;
    LS_s = '--';
    anchor_vis = 'off';
    LW_r = 2;
    LS_r = '-';
  end
  
  pos_tmp = h.data.pos_s*3*h.data.ymax/h.paras.ptlength;
  set(h.pos,'Ydata',pos_tmp,'LineStyle',LS_s,'LineWidth',LW_s)
  
  pos_r_tmp = h.data.pos_r*3*h.data.ymax/h.paras.ptlength;
  set(h.pos_r,'Ydata',pos_r_tmp,'LineStyle',LS_r,'LineWidth',LW_r)
  
  set(h.scatter_pos0,'visible',anchor_vis)
  set(h.scatter_speed0,'visible',anchor_vis)
  
  set(h.scatter_pos,'Xdata',h.data.bh.time(h.data.idx_align_pos(h.data.idx_align_pos_mask)+1))
  set(h.scatter_pos,'Ydata',ones(1,sum(h.data.idx_align_pos_mask)),'visible',anchor_vis);
  
  set(h.scatter_speed,'Xdata',h.data.bh.time(h.data.idx_align_speed(h.data.idx_align_speed_mask)+1))
  set(h.scatter_speed,'Ydata',-ones(1,sum(h.data.idx_align_speed_mask)),'visible',anchor_vis);
  
  
  
  err = h.data.pos_s(1:h.data.bh.datapoints)-(h.data.bh.pos(1:h.data.bh.datapoints) + h.paras.pos_offset);
  err = min([abs(err);abs(err-h.paras.pos_offset)])/(h.paras.pos_offset/100);
  
  err_r = h.data.pos_r(1:h.data.bh.datapoints)-(h.data.bh.pos(1:h.data.bh.datapoints) + h.paras.pos_offset);
  err_r = min([abs(err_r);abs(err_r-h.paras.pos_offset)])/(h.paras.pos_offset/100);
  
  
  
  set(h.err,'Ydata',err,'LineStyle',LS_s,'LineWidth',LW_s)
  set(h.err_r,'Ydata',err_r,'LineStyle',LS_r,'LineWidth',LW_r)
  
  err_ylim = max(10,max(imgaussfilt(err,10)));
  ylim(h.ax3,[-err_ylim,err_ylim])
  
  if get(h.radio_stops,'Value')
    cum_pos = zeros(size(h.data.bh.time));
    cum_pos(h.data.idx_align_pos(h.data.idx_align_pos_mask)+1) = 1;
    cum_pos = cumsum(cum_pos);
    
    cum_speed = zeros(size(h.data.bh.time));
    cum_speed(h.data.idx_align_speed(h.data.idx_align_speed_mask)+1) = 1;
    cum_speed = cumsum(cum_speed);
    
    cum_diff = cum_pos-cum_speed;
    
    ylim(h.ax2,[-5,5])
    
    set(h.ax2_err,'Xdata',h.data.bh.time)
    set(h.ax2_err,'Ydata',cum_diff)
    set(get(h.ax2,'YLabel'),'String','diff. anchor points')
  else
    err_r_0 = h.data.pos_r(1:h.data.bh.datapoints)-h.data.pos_s(1:h.data.bh.datapoints);
    err_r_0 = min([abs(err_r_0);abs(err_r_0-h.paras.pos_offset)])/(h.paras.pos_offset/100);
    
    err_ylim = max(10,max(imgaussfilt(err_r_0,10)));
    ylim(h.ax2,[-err_ylim,err_ylim])
    
    set(h.ax2_err,'Xdata',h.data.time)
    set(h.ax2_err,'Ydata',err_r_0)
    set(get(h.ax2,'YLabel'),'String','TTL vs stop')
  end
    
  
% --- Outputs from this function are returned to the command line.
function varargout = alignData_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = h.output;
show_info();


function remove_align_point(hObject,eventdata,ax)
  
  pos = get(ax,'CurrentPoint');
  
  x = pos(1,1);
  y = round(pos(1,2));
  
  h = guidata(hObject);
  
  if y>0 && y<h.data.ymax/3
    [min_val,idx_rm] = min(abs(h.data.bh.time(h.data.idx_align_pos+1)-x));
    if min_val < 5
      h.data.idx_align_pos_mask(idx_rm) = ~h.data.idx_align_pos_mask(idx_rm);
      h = run_align(h);
      guidata(hObject,h);
      
      update_plot(h)
    end
  elseif y<0 && y>-h.data.ymax/3
    [min_val,idx_rm] = min(abs(h.data.bh.time(h.data.idx_align_speed+1)-x));
    if min_val < 5
      h.data.idx_align_speed_mask(idx_rm) = ~h.data.idx_align_speed_mask(idx_rm);
      h = run_align(h);
      guidata(hObject,h);
      update_plot(h)
    end
  end
  
  
  
% --- Executes on button press in button_reset.
function button_reset_Callback(hObject, eventdata, h)
% hObject    handle to button_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  h = set_paras(h);
  if h.nSessions
    bhfile=dir(pathcat(h.pathIn,h.folders(h.currSession).name,'*m.txt'));
    h = read_bh_file(h,pathcat(h.pathIn,h.folders(h.currSession).name,bhfile.name));
  else
%      h.pathIn
%      bhfile.name
%      bhfile=dir(pathcat(h.pathIn,'*m.txt'));
    h = read_bh_file(h,h.pathIn);
  end
  
  h = run_methods(h);
  guidata(hObject,h);
  
  

% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, h)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  save_data(h)
  
  h = updateSaveStatus(h,true);
  guidata(hObject,h);
  
  
  
function save_data(h)
  
  alignedData = struct;
  
  if get(h.radio_stops,'Value')
    alignedData.position = h.data.pos_s(1:h.data.bh.datapoints);            % position
  else
    alignedData.position = h.data.pos_r(1:h.data.bh.datapoints);
  end
  alignedData.speed = h.data.speed;
  alignedData.time = h.data.time;        % time
  alignedData.frame = h.data.frame_num;
  
  alignedData.resampled = struct;
  
  %% resample stuff at 15 Hz by averaging
  for i=1:h.paras.nframe
    idx = find(h.data.frame_num==i);
    alignedData.resampled.time(i) = mean(h.data.time(idx));        % time
    alignedData.resampled.frame(i) = i;                            % frame number
    
    if get(h.radio_stops,'Value')
      val_tmp = sort(h.data.pos_s(idx));
      alignedData.resampled.position(i) = val_tmp(round(length(val_tmp)/2));   % quasi median filtering (avoiding half values!)
    else
      val_tmp = sort(h.data.pos_r(idx));
      alignedData.resampled.position(i) = val_tmp(round(length(val_tmp)/2));
    end
    
    alignedData.resampled.speed(i) = mean(h.data.speed(idx));      % speed
  end
  alignedData.resampled.duration = alignedData.resampled.time(end);
  
  h.paras.binwidth = h.paras.ptlength/h.paras.nbin;
  alignedData.resampled.binpos = min(floor(alignedData.resampled.position/h.paras.binwidth)+1,h.paras.nbin);
    
  %% defining run epochs by finding super-threshold movement speed
  sm = ones(1,10+1);
  alignedData.resampled.runrest = imerode(imdilate(alignedData.resampled.speed >= h.paras.runthres,sm),sm);   %% filling holes of size up to 5 (=1/3 sec)
  alignedData.resampled.longrunperiod = bwareaopen(alignedData.resampled.runrest,h.paras.lr_min);    % find lonrunperiods as connected areas > lr_min of size
  
  %% create dwell time histogram
  alignedData.resampled.dwelltime = zeros(1,h.paras.nbin);
  binpos = max(1,alignedData.resampled.binpos(alignedData.resampled.longrunperiod));
  for i=1:length(binpos)
    alignedData.resampled.dwelltime(binpos(i)) = alignedData.resampled.dwelltime(binpos(i)) + 1/h.paras.f;
  end
  alignedData.resampled.norm_dwelltime = alignedData.resampled.dwelltime/sum(alignedData.resampled.dwelltime);
  
  [~,fileSave,~] = fileparts(h.pathBHfile);
  
  pathSave = pathcat(h.pathSession,sprintf('%s_aligned.mat',fileSave));
  save(pathSave,'alignedData','-v7.3')
  
  if get(h.checkbox_savefigure,'value')
    pathSave = pathcat(h.pathSession,sprintf('%s_aligned.png',fileSave));
    export_fig(gcf,pathSave,'-png')
    disp(sprintf('figure saved as %s',pathSave))
  end
  

function h = updateSaveStatus(h,new_status)
  
  h.save_status = new_status;
  if h.save_status
    set(h.button_save,'BackgroundColor','green');
  else
    set(h.button_save,'BackgroundColor','red');
  end
  


% --- Executes on button press in button_next.
function button_next_Callback(hObject, eventdata, h)
% hObject    handle to button_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  if ~h.nSessions
    
    if get(h.checkbox_autosave,'Value')
      save_data(h)
    end
    
    [fileName, h.pathSession] = uigetfile({'*.txt'},'Choose behavior file',h.pathSession);
    h.pathIn = pathcat(h.pathSession,fileName);
    h.pathBHfile = h.pathIn;
    h = read_bh_file(h,h.pathBHfile);
    h = run_methods(h);
    h = updateSaveStatus(h,false);
    
    guidata(hObject,h);
    
  elseif h.currSession < h.nSessions
    if get(h.checkbox_autosave,'Value')
      save_data(h)
    end
    h.save_status = false;
    h.currSession = h.currSession + 1;
    if h.currSession == h.nSessions
      set(h.button_next,'string','Finish')
    else
      set(h.button_previous,'enable','on')
    end
    h.pathSession = pathcat(h.pathIn,h.folders(h.currSession).name);
    bhfile=dir(pathcat(h.pathIn,h.folders(h.currSession).name,'*m.txt'));
    h.pathBHfile = pathcat(h.pathIn,h.folders(h.currSession).name,bhfile.name);
    h = read_bh_file(h,h.pathBHfile);
%      h = read_bh_file(h,pathcat(h.pathIn,h.folders(h.currSession).name,bhfile.name));
    h = run_methods(h);
    guidata(hObject,h);
  elseif h.currSession == h.nSessions
    if get(h.checkbox_autosave,'Value')
      save_data(h)
    end
    delete(gcf)
  end
    


% --- Executes on button press in button_previous.
function button_previous_Callback(hObject, eventdata, h)
% hObject    handle to button_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  if h.currSession > 1
    if get(h.checkbox_autosave,'Value')
      save_data(h)
    end
    h.save_status = false;
    h.currSession = h.currSession - 1;
    if h.currSession == 1
      set(h.button_previous,'enable','off')
    end
    h.pathSession = pathcat(h.pathIn,h.folders(h.currSession).name);
    bhfile=dir(pathcat(h.pathIn,h.folders(h.currSession).name,'*m.txt'));
    h.pathBHfile = pathcat(h.pathIn,h.folders(h.currSession).name,bhfile.name);
    h = read_bh_file(h,h.pathBHfile);
    h = run_methods(h);
    guidata(hObject,h);
  end
  
  
  
  

function entry_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to entry_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_smooth as text
%        str2double(get(hObject,'String')) returns contents of entry_smooth as a double


% --- Executes during object creation, after setting all properties.
function entry_smooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entry_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close alignData.
function alignData_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to alignData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);



function entry_tracklength_Callback(hObject, eventdata, handles)
% hObject    handle to entry_tracklength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_tracklength as text
%        str2double(get(hObject,'String')) returns contents of entry_tracklength as a double


% --- Executes during object creation, after setting all properties.
function entry_tracklength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entry_tracklength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entry_pos_offset_Callback(hObject, eventdata, handles)
% hObject    handle to entry_pos_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_pos_offset as text
%        str2double(get(hObject,'String')) returns contents of entry_pos_offset as a double


% --- Executes during object creation, after setting all properties.
function entry_pos_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entry_pos_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function entry_col_reward_Callback(hObject, eventdata, handles)
% hObject    handle to entry_col_reward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_col_reward as text
%        str2double(get(hObject,'String')) returns contents of entry_col_reward as a double


% --- Executes during object creation, after setting all properties.
function entry_col_reward_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entry_col_reward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entry_col_TTL_Callback(hObject, eventdata, handles)
% hObject    handle to entry_col_TTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_col_TTL as text
%        str2double(get(hObject,'String')) returns contents of entry_col_TTL as a double


% --- Executes during object creation, after setting all properties.
function entry_col_TTL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entry_col_TTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entry_reward_position_Callback(hObject, eventdata, handles)
% hObject    handle to entry_reward_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entry_reward_position as text
%        str2double(get(hObject,'String')) returns contents of entry_reward_position as a double


% --- Executes during object creation, after setting all properties.
function entry_reward_position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entry_reward_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_autosave.
function checkbox_autosave_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_autosave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_autosave


% --- Executes on button press in radio_TTL.
function radio_TTL_Callback(hObject, eventdata, h)
% hObject    handle to radio_TTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_TTL
   
  set(h.radio_stops,'Value',~get(h.radio_stops,'Value'))
  update_plot(h)
  
% --- Executes on button press in radio_stops.
function radio_stops_Callback(hObject, eventdata, h)
% hObject    handle to radio_stops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_stops
  
  set(h.radio_TTL,'Value',~get(h.radio_TTL,'Value'))
  update_plot(h)
