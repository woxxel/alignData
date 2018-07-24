%% GUI to align position and speed data
%%
%% INPUT via getappdata from "prep_behaviour.m": preloaded data from behavior-text-file
%% OUTPUT via setappdata to "prep_behavior.m": aligned position values

%% written 24.07.2018 by A. Schmidt
%%
%% TODOs:
%% - make loading possible from GUI itself
%% - make program independent from other methods
%% - give estimate, when data seems to be alright (from overall number of alignment points?
%% - make more descriptive?
%% - automatically find number of columns in text-file



function varargout = pos_speed_align(varargin)
% POS_SPEED_ALIGN MATLAB code for pos_speed_align.fig
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
%      applied to the GUI before pos_speed_align_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pos_speed_align_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pos_speed_align

% Last Modified by GUIDE v2.5 23-Jul-2018 18:55:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pos_speed_align_OpeningFcn, ...
                   'gui_OutputFcn',  @pos_speed_align_OutputFcn, ...
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


% --- Executes just before pos_speed_align is made visible.
function pos_speed_align_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pos_speed_align (see VARARGIN)

  handles.status = 0;
  
  handles.data = struct;
  handles.data.bh = getappdata(0,'behavior_data');
  handles.paras = getappdata(0,'parameter');
  
  handles.data.speed = handles.data.bh.speed(1:handles.data.bh.datapoints);
  handles.data.time = handles.data.bh.time(1:handles.data.bh.datapoints);
  handles.data.duration = handles.data.time(end);
  
  handles = smooth_data(handles);
  guidata(hObject,handles);
  
  handles = run_align(handles);
  guidata(hObject,handles);
  
  handles = plot_align(handles);
  guidata(hObject,handles);
  
  
  % Choose default command line output for pos_speed_align
  handles.output = hObject;

  % Update handles structure
  guidata(hObject, handles);

% UIWAIT makes pos_speed_align wait for user response (see UIRESUME)
% uiwait(handles.alignData);


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
    
  
  
function h = run_align(h)
  
  h.data.pos = zeros(1,length(h.data.bh.pos));
  offset = 1;
  i = 2;
  
  idx_align_pos = h.data.idx_align_pos(h.data.idx_align_pos_mask);
  idx_align_speed = h.data.idx_align_speed(h.data.idx_align_speed_mask);
  
  while offset < h.data.bh.datapoints
    n_loc = idx_align_pos(i)-idx_align_pos(i-1);
    n_speed = idx_align_speed(i)-idx_align_speed(i-1);
    %% check if distances of stopping periods are the same
    if n_loc~=n_speed   %% if not, resample location points and write to new position
      loc_tmp = h.data.bh.pos(idx_align_pos(i-1)+1:idx_align_pos(i));
      
      time_speed = linspace(h.data.bh.time(idx_align_speed(i-1)+1),h.data.bh.time(idx_align_speed(i)),n_speed);
      time_loc = linspace(h.data.bh.time(idx_align_speed(i-1)+1),h.data.bh.time(idx_align_speed(i)),n_loc);
      
      h.data.pos(offset:offset+n_speed-1) = interp1(time_loc,loc_tmp,time_speed);
      
      offset = offset + n_speed;
    else      %% if they are, just write old to new data
      loc_tmp = h.data.bh.pos(idx_align_pos(i-1)+1:idx_align_pos(i));
      h.data.pos(offset:offset+n_speed-1) = loc_tmp;
      
      offset = offset + length(loc_tmp);
    end
    
    %% if there is no stopping period after recording stops, assign remaining positions linearly
    if i == length(idx_align_pos) || i == length(idx_align_speed)
      n_speed = min(h.data.bh.datapoints - idx_align_speed(i),h.data.bh.datapoints - idx_align_pos(i));
      
      if n_speed > 0
        h.data.pos(offset:offset+n_speed-1) = h.data.bh.pos(idx_align_pos(i)+1:idx_align_pos(i)+n_speed);
      end
      break
    end
    i = i+1;
  end
  
  
  
  
function h = plot_align(h)
  
  cla(h.ax,'reset')
  cla(h.ax2,'reset')
  cla(h.ax3,'reset')
  
  h.data.ymax = ceil(max(h.data.speed));
  
  plot(h.ax,h.data.time,-imgaussfilt(h.data.speed,10),'b','Hittest','off')
  
  hold(h.ax,'on')
  bar(h.ax,h.data.bh.time,~h.data.speed_bool*3*h.data.ymax,1,'FaceColor',[0.8 0.8 0.8],'Hittest','off')
  
  pos_tmp = h.data.pos*3*h.data.ymax/h.paras.totallength;
  h.pos = plot(h.ax,h.data.bh.time,pos_tmp,'k','Hittest','off','DisplayName','animal position (aligned)');
  h.pos0 = plot(h.ax,h.data.time,h.data.bh.pos(1:h.data.bh.datapoints)*3*h.data.ymax/h.paras.totallength,'r--','Hittest','off','DisplayName','animal position');
  
  h.reward = plot(h.ax,[0,h.data.duration],[60,60]*3*h.data.ymax/h.paras.nbin,'g--','LineWidth',2,'Hittest','off','DisplayName','reward position');
  h.gate = plot(h.ax,[0,h.data.duration],[20,20]*3*h.data.ymax/h.paras.nbin,'r--','LineWidth',2,'Hittest','off','DisplayName','gate position');
  
  plot(h.ax,[h.data.duration, h.data.duration],[-h.data.ymax,3*h.data.ymax],'k--')
  
  scatter(h.ax,h.data.bh.time(h.data.idx_align_pos+1),ones(1,length(h.data.idx_align_pos)),'ko','Hittest','off')
  scatter(h.ax,h.data.bh.time(h.data.idx_align_speed+1),-ones(1,length(h.data.idx_align_speed)),'ro','Hittest','off')
  
  h.scatter_pos = scatter(h.ax,h.data.bh.time(h.data.idx_align_pos(h.data.idx_align_pos_mask)+1),ones(1,sum(h.data.idx_align_pos_mask)),'ko','filled','Hittest','off','DisplayName','stops (position)');
  h.scatter_speed = scatter(h.ax,h.data.bh.time(h.data.idx_align_speed(h.data.idx_align_speed_mask)+1),-ones(1,sum(h.data.idx_align_speed_mask)),'ro','filled','Hittest','off','DisplayName','stops (speed)');
  
  hold(h.ax,'off')
  ylim(h.ax,[-h.data.ymax 3*h.data.ymax])
  yticks(h.ax,[-h.data.ymax,0,h.data.ymax])
  yticklabels(h.ax,[h.data.ymax,0,1])
  ylabel(h.ax,'Position')
  xlim(h.ax,[0,h.data.bh.time(end)])
  
  
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
  h.diff = plot(h.ax2,h.data.bh.time,cum_diff,'r');
  plot(h.ax2,[h.data.duration, h.data.duration],[-5,5],'k--')
  hold(h.ax2,'off')
  xlim(h.ax2,[0,h.data.bh.time(end)])
  ylim(h.ax2,[-5,5])
  ylabel(h.ax2,'diff. anchor points')
  
  err = h.data.pos(1:h.data.bh.datapoints)-h.data.bh.pos(1:h.data.bh.datapoints);
  err = min([abs(err);abs(err-h.paras.pos_offset)])/(h.paras.pos_offset/100);
%    err
  hold(h.ax3,'on')
  h.err = plot(h.ax3,h.data.time,err,'r');
  
  err_ylim = max(10,max(imgaussfilt(err,10)));
  plot(h.ax3,[h.data.duration, h.data.duration],[-err_ylim,err_ylim],'k--')
  hold(h.ax3,'off')
  xlim(h.ax3,[0,h.data.bh.time(end)])
  
  ylim(h.ax3,[-err_ylim,err_ylim])
  xlabel(h.ax3,'time [s]')
  ylabel(h.ax3,'deviation [cm]')
  
%    plot cumulative points or rather difference of cummulative points (to see, where to remove etc)
  
  
  linkaxes([h.ax,h.ax2,h.ax3],'x')
  
%    linkdata on
  
  
function update_plot(h)
  
  pos_tmp = h.data.pos*3*h.data.ymax/h.paras.totallength;
  set(h.pos,'Ydata',pos_tmp)
  
  set(h.scatter_pos,'Xdata',h.data.bh.time(h.data.idx_align_pos(h.data.idx_align_pos_mask)+1))
  set(h.scatter_pos,'Ydata',ones(1,sum(h.data.idx_align_pos_mask)));
  
  set(h.scatter_speed,'Xdata',h.data.bh.time(h.data.idx_align_speed(h.data.idx_align_speed_mask)+1))
  set(h.scatter_speed,'Ydata',-ones(1,sum(h.data.idx_align_speed_mask)));
  
  cum_pos = zeros(size(h.data.bh.time));
  cum_pos(h.data.idx_align_pos(h.data.idx_align_pos_mask)+1) = 1;
  cum_pos = cumsum(cum_pos);
  
  cum_speed = zeros(size(h.data.bh.time));
  cum_speed(h.data.idx_align_speed(h.data.idx_align_speed_mask)+1) = 1;
  cum_speed = cumsum(cum_speed);
  
  cum_diff = cum_pos-cum_speed;
  
  set(h.diff,'Ydata',cum_diff)
  
  err = h.data.pos(1:h.data.bh.datapoints)-h.data.bh.pos(1:h.data.bh.datapoints);
  err = min([abs(err);abs(err-h.paras.pos_offset)])/(h.paras.pos_offset/100);
  
  set(h.err,'Ydata',err)
  err_ylim = max(10,max(imgaussfilt(err,10)));
  ylim(h.ax3,[-err_ylim,err_ylim])
%    refreshdata(h.pos)
%    refreshdata(h.scatter_pos)
%    refreshdata(h.scatter_speed)
  
  
  
  
  
  

% --- Outputs from this function are returned to the command line.
function varargout = pos_speed_align_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
function button_reset_Callback(hObject, eventdata, handles)
% hObject    handle to button_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  handles = smooth_data(handles);
  guidata(hObject,handles);
  
  handles = run_align(handles);
  guidata(hObject,handles);
  
  handles = plot_align(handles);
  guidata(hObject,handles);
  
  

% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
  setappdata(0,'new_position',handles.data.pos)
  
  delete(gcf)

% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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
