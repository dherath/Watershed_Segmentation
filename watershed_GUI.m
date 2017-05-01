function varargout = watershed_GUI(varargin)
clc;
% WATERSHED_GUI MATLAB code for watershed_GUI.fig
%      WATERSHED_GUI, by itself, creates a new WATERSHED_GUI or raises the existing
%      singleton*.
%
%      H = WATERSHED_GUI returns the handle to a new WATERSHED_GUI or the handle to
%      the existing singleton*.
%
%      WATERSHED_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WATERSHED_GUI.M with the given input arguments.
%
%      WATERSHED_GUI('Property','Value',...) creates a new WATERSHED_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before watershed_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to watershed_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help watershed_GUI

% Last Modified by GUIDE v2.5 10-Jan-2016 06:21:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @watershed_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @watershed_GUI_OutputFcn, ...
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


% --- Executes just before watershed_GUI is made visible.
function watershed_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to watershed_GUI (see VARARGIN)

% Choose default command line output for watershed_GUI
handles.output = hObject;

% Update handles structure
imshow(zeros(256,256));axis off;axis image;
guidata(hObject, handles);

% UIWAIT makes watershed_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = watershed_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('rawdat.mat');
I=uint8(rawdat);
[x,y,z]=size(I);

%size variables---------------------------------
    handles.I=I;
    handles.final=I;
    handles.x=x;
    handles.y=y;
    handles.z=z;
    %------defaut view is from side------------
    handles.sv=1;%side view->1
    handles.fv=0;%front view->0
    handles.tv=0;%top view->0
    segment=0;
    handles.segment=double(segment);
    final3D=uint8(zeros(x,y,z));
    handles.Lrgb_sv=final3D;
    handles.Lrgb_tv=final3D;%uint8(zeros(256,180,256));
    handles.Lrgb_fv=final3D;%uint8(zeros(256,180,256));
    guidata(hObject,handles);
    %------------------------------------------
%     %handling important true/false
%     flag.segment=0;
%     guidata(hObject,flag);
%--------------------------------
k=I(:,:,1);
axes(handles.axes1);
imshow(k);

function segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function manual_Callback(hObject, eventdata, handles)
% hObject    handle to manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.final=handles.I;
p=round(handles.slider_val);%position
I=uint8(handles.I);

countk=0;
count1=0;
count2=[0,0,0];
%-------based on viewpoint--------------------------------------------
if handles.sv==1 
    k=I(:,:,p);
    limit=180;
elseif handles.fv==1
    k=squeeze(I(p,:,:));
    limit=256;
elseif handles.tv==1
    k=squeeze(I(:,p,:));
    limit=256;
end
axes(handles.axes1);
temp=uint8(roipoly(k));
%---------getting mask---------------------------
mask=logical(imdilate(temp,strel('disk',5)));
%---------finding mean-intensity region----------------
   temp=temp.*k;
   pixels=nnz(temp);
   area=double(sum(sum(temp))/(pixels));
%--------getting region-------------------------
xs=handles.x;
ys=handles.y;
zs=handles.z;
final3D=uint8(zeros(xs,ys,zs));
Lrgb_sv=final3D;
Lrgb_tv=final3D;%uint8(zeros(256,180,256));
Lrgb_fv=final3D;%uint8(zeros(256,180,256));
for i=1:xs
    for j=1:ys
        for cn=1:zs
            if I(i,j,cn)>area*0.5 && I(i,j,cn)<area*1.5
                final3D(i,j,cn)=1;
            end
        end
    end   
end
disp('final pixel mask');
disp(nnz(final3D));
%%-------------------------------------
while nnz(mask)>pixels*0.5
mask=imerode(mask,strel('disk',5));
end
disp('global point');disp(nnz(mask));
global_point=mask;%propogated to next slide
for count=1:3
   range=imdilate(global_point,strel('disk',25)); 
end
disp('dilated range mask');disp(nnz(range));
for no=1:limit
    countk=countk+1;
    if handles.sv==1
        
        final=squeeze(final3D(:,:,no));
        mask=imreconstruct(global_point,logical(final));
        k=uint8(squeeze(I(:,:,no)));
        final=uint8(mask).*uint8(range);
       final=imerode(final,[1,1;1,1]);
       final=bwareaopen(logical(final),round(pixels*0.3));
       final=uint8(final).*uint8(k);
       final3D(:,:,no)=final;

    elseif handles.tv==1
        
        final=squeeze(final3D(:,no,:));
        mask=imreconstruct(global_point,logical(final));
        k=uint8(squeeze(I(:,no,:)));
        final=uint8(mask).*uint8(range);
       final=imerode(final,[1,1;1,1]);
        final=bwareaopen(logical(final),round(pixels*0.3));
        final=uint8(final).*uint8(k);
        final3D(:,no,:)=final;
    elseif handles.fv==1
        
        final=squeeze(final3D(no,:,:));
        mask=imreconstruct(global_point,logical(final));
        k=uint8(squeeze(I(no,:,:)));
        final=uint8(mask).*uint8(range);
       final=imerode(final,[1,1;1,1]);
       final=bwareaopen(logical(final),round(pixels*0.3));
        final=uint8(final).*uint8(k);
        final3D(no,:,:)=final;
    end   
end
disp('Actual mask');
disp(nnz(final3D));
%------------watershed-----------------
hy=fspecial('sobel');
hx=hy';
for t=1:3
    count1=count1+1;
    if t==1
        slices=180;
    elseif t==2
        slices=256;
    elseif t==3
        slices=256;
    end
    for tot=1:slices
        count2(1,t)=count2(1,t)+1;
         if t==1
            final=squeeze(final3D(:,:,tot));
         elseif t==2        
            final=squeeze(final3D(:,tot,:));
         elseif t==3
             final=squeeze(final3D(tot,:,:));
         end
        Gy=imfilter(double(final),hy,'replicate');
        Gx=imfilter(double(final),hx,'replicate');
        gradmag=sqrt(Gx.^2+Gy.^2);
        
        
        se=strel('disk',5);
        Ie=imerode(final,se);
        Iobr=imreconstruct(Ie,final);
        
        Iobrd=imdilate(Iobr,se);
        Iobrcbr=imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
        Iobrcbr=imcomplement(Iobrcbr);

        fgm=imregionalmax(Iobrcbr);
        I2=final;
        I2(fgm)=255;

        se2=strel(ones(3,3));
        fgm2=imclose(fgm,se2);
        fgm3=imerode(fgm2,se2);
        fgm4=bwareaopen(fgm3,50);
        I3=Iobrcbr;
        I3(fgm4)=255;

        bw=im2bw(Iobrcbr(:,:,1),graythresh(Iobrcbr(:,:,1)));

        D=bwdist(bw);
        
        ors=D|fgm4;
        ors=double(ors);
        gradmag2=imimposemin(gradmag,ors);
        
        L=watershed(gradmag2);
        %Lrgb=label2rgb(L,'jet','w','shuffle');

         if t==1
           Lrgb_sv(:,:,tot)=L;%side view->180 slices
         elseif t==2        
           Lrgb_tv(:,tot,:)=L;%top view->256
         elseif t==3
           Lrgb_fv(tot,:,:)=L;%front view->256
         end      
    end    
end
L_tot=(Lrgb_sv+Lrgb_tv+Lrgb_fv)./3;%new line added
handles.Lrgb_sv=Lrgb_sv;
handles.Lrgb_fv=Lrgb_fv;
handles.Lrgb_tv=Lrgb_tv;
handles.L_tot=L_tot;
segment=100;
handles.segment=double(segment);
guidata(hObject,handles);
disp('process worked');
if handles.sv==1
	img=label2rgb(squeeze(Lrgb_sv(:,:,p)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(I(:,:,p)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);
elseif handles.tv==1
    img=label2rgb(squeeze(Lrgb_tv(:,p,:)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(I(:,p,:)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);
elseif handles.fv==1
    img=label2rgb(squeeze(Lrgb_fv(p,:,:)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(I(p,:,:)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);   
end
disp(handles.segment);

% --------------------------------------------------------------------
function automatic_Callback(hObject, eventdata, handles)
% hObject    handle to automatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.final=handles.I;
p=round(handles.slider_val);%position
I=uint8(handles.I);
%----------predefined masks------------------------------------------
load('b_sv.mat');
load('b_tv.mat');
load('b_fv.mat');
%--------------------------------------------------------------------
countk=0;
count1=0;
count2=[0,0,0];
%-------based on viewpoint--------------------------------------------
if handles.sv==1 
    k=I(:,:,p);
    limit=180;
    temp=b_sv(:,:,p);
elseif handles.fv==1
    k=squeeze(I(p,:,:));
    limit=256;
    temp=squeeze(b_fv(p,:,:));
elseif handles.tv==1
    k=squeeze(I(:,p,:));
    limit=256;
    temp=squeeze(b_tv(:,p,:));
end
%axes(handles.axes1);
%temp=uint8(roipoly(k));
disp(size(temp));
disp(size(k));
%---------getting mask---------------------------
mask=logical(imdilate(temp,strel('disk',5)));
%---------finding mean-intensity region----------------
   temp=temp.*k;
   pixels=nnz(temp);
   area=double(sum(sum(temp))/(pixels));
%--------getting region-------------------------
xs=handles.x;
ys=handles.y;
zs=handles.z;
final3D=uint8(zeros(xs,ys,zs));
Lrgb_sv=final3D;
Lrgb_tv=final3D;%uint8(zeros(256,180,256));
Lrgb_fv=final3D;%uint8(zeros(256,180,256));
for i=1:xs
    for j=1:ys
        for cn=1:zs
            if I(i,j,cn)>area*0.5 && I(i,j,cn)<area*1.5
                final3D(i,j,cn)=1;
            end
        end
    end   
end
disp('final pixel mask');
disp(nnz(final3D));
%%-------------------------------------
while nnz(mask)>pixels*0.5
mask=imerode(mask,strel('disk',5));
end
disp('global point');disp(nnz(mask));
global_point=mask;%propogated to next slide
for count=1:3
   range=imdilate(global_point,strel('disk',25)); 
end
disp('dilated range mask');disp(nnz(range));
for no=1:limit
    countk=countk+1;
    if handles.sv==1
        
        final=squeeze(final3D(:,:,no));
        mask=imreconstruct(global_point,logical(final));
        k=uint8(squeeze(I(:,:,no)));
        final=uint8(mask).*uint8(range);
       final=imerode(final,[1,1;1,1]);
       final=bwareaopen(logical(final),round(pixels*0.3));
       final=uint8(final).*uint8(k);
       final3D(:,:,no)=final;

    elseif handles.tv==1
        
        final=squeeze(final3D(:,no,:));
        mask=imreconstruct(global_point,logical(final));
        k=uint8(squeeze(I(:,no,:)));
        final=uint8(mask).*uint8(range);
       final=imerode(final,[1,1;1,1]);
        final=bwareaopen(logical(final),round(pixels*0.3));
        final=uint8(final).*uint8(k);
        final3D(:,no,:)=final;
    elseif handles.fv==1
        
        final=squeeze(final3D(no,:,:));
        mask=imreconstruct(global_point,logical(final));
        k=uint8(squeeze(I(no,:,:)));
        final=uint8(mask).*uint8(range);
       final=imerode(final,[1,1;1,1]);
       final=bwareaopen(logical(final),round(pixels*0.3));
        final=uint8(final).*uint8(k);
        final3D(no,:,:)=final;
    end   
end
disp('Actual mask');
disp(nnz(final3D));
%------------watershed-----------------
hy=fspecial('sobel');
hx=hy';
for t=1:3
    count1=count1+1;
    if t==1
        slices=180;
    elseif t==2
        slices=256;
    elseif t==3
        slices=256;
    end
    for tot=1:slices
        count2(1,t)=count2(1,t)+1;
         if t==1
            final=squeeze(final3D(:,:,tot));
         elseif t==2        
            final=squeeze(final3D(:,tot,:));
         elseif t==3
             final=squeeze(final3D(tot,:,:));
         end
        Gy=imfilter(double(final),hy,'replicate');
        Gx=imfilter(double(final),hx,'replicate');
        gradmag=sqrt(Gx.^2+Gy.^2);
        
        
        se=strel('disk',5);
        Ie=imerode(final,se);
        Iobr=imreconstruct(Ie,final);
        
        Iobrd=imdilate(Iobr,se);
        Iobrcbr=imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
        Iobrcbr=imcomplement(Iobrcbr);

        fgm=imregionalmax(Iobrcbr);
        I2=final;
        I2(fgm)=255;

        se2=strel(ones(3,3));
        fgm2=imclose(fgm,se2);
        fgm3=imerode(fgm2,se2);
        fgm4=bwareaopen(fgm3,50);
        I3=Iobrcbr;
        I3(fgm4)=255;

        bw=im2bw(Iobrcbr(:,:,1),graythresh(Iobrcbr(:,:,1)));

        D=bwdist(bw);
        
        ors=D|fgm4;
        ors=double(ors);
        gradmag2=imimposemin(gradmag,ors);
        
        L=watershed(gradmag2);
        %Lrgb=label2rgb(L,'jet','w','shuffle');

         if t==1
           Lrgb_sv(:,:,tot)=L;%side view->180 slices
         elseif t==2        
           Lrgb_tv(:,tot,:)=L;%top view->256
         elseif t==3
           Lrgb_fv(tot,:,:)=L;%front view->256
         end      
    end    
end
L_tot=(Lrgb_sv+Lrgb_tv+Lrgb_fv)./3;%new line added
handles.Lrgb_sv=Lrgb_sv;
handles.Lrgb_fv=Lrgb_fv;
handles.Lrgb_tv=Lrgb_tv;
handles.L_tot=L_tot;
segment=100;
handles.segment=double(segment);
guidata(hObject,handles);
disp('process worked');
if handles.sv==1
	img=label2rgb(squeeze(Lrgb_sv(:,:,p)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(I(:,:,p)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);
elseif handles.tv==1
    img=label2rgb(squeeze(Lrgb_tv(:,p,:)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(I(:,p,:)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);
elseif handles.fv==1
    img=label2rgb(squeeze(Lrgb_fv(p,:,:)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(I(p,:,:)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);   
end
disp(handles.segment);


% --- Executes on slider movement.

function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
disp('segment value->');
disp(handles.segment);
cla;
if handles.sv==1
     if hObject.Value > 180
     hObject.Value=handles.slider_val;    
     end
    hObject.Max=handles.z;
    val=get(hObject,'Value');
    if handles.segment==0
        axes(handles.axes1);
        imshow(handles.final(:,:,round(val)));
    elseif handles.segment==100
        img=label2rgb(squeeze(handles.Lrgb_sv(:,:,round(val))));
        axes(handles.axes1);
        imshow(squeeze(handles.I(:,:,round(val))));axis off;axis image;
        hold on;
        himage=imshow(img);axis off;axis image;
        set(himage,'Alphadata',0.3);
    end
elseif handles.fv==1
    hObject.Max=handles.y;
    val=get(hObject,'Value');
    if handles.segment==0  
        axes(handles.axes1);
        imshow(squeeze(handles.final(round(val),:,:)));
    elseif handles.segment==100
        img=label2rgb(squeeze(handles.Lrgb_fv(round(val),:,:)));
        axes(handles.axes1);
        imshow(squeeze(handles.I(round(val),:,:)));axis off;axis image
        hold on;
        himage=imshow(img);axis off;axis image;
        set(himage,'Alphadata',0.3);
    end
elseif handles.tv==1
    hObject.Max=handles.x;
    val=get(hObject,'Value');
    if handles.segment==0
        axes(handles.axes1);
        imshow(squeeze(handles.final(:,round(val),:)));
    elseif handles.segment==100
        img=label2rgb(squeeze(handles.Lrgb_tv(:,round(val),:)));
        axes(handles.axes1);
        imshow(squeeze(handles.I(:,round(val),:)));axis off;axis image;
        hold on;
        himage=imshow(img);axis off;axis image;
        set(himage,'Alphadata',0.3);
    end
end
handles.limit=hObject.Max;
handles.slider_val=get(hObject,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function viewpoint_Callback(hObject, eventdata, handles)
% hObject    handle to viewpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function side_view_Callback(hObject, eventdata, handles)
% hObject    handle to side_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fv=0;
handles.tv=0;
handles.sv=1;
handles.slider_val=50;
handles.xs=handles.x;
handles.ys=handles.y;
guidata(hObject,handles);
    cla;
if handles.segment==100
	img=label2rgb(squeeze(handles.Lrgb_sv(:,:,50)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(handles.I(:,:,50)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);
else
    axes(handles.axes1);
    imshow(squeeze(handles.final(:,:,1)));axis off;axis image;
end

% --------------------------------------------------------------------
function front_view_Callback(hObject, eventdata, handles)
% hObject    handle to front_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fv=1;
handles.tv=0;
handles.sv=0;
handles.slider_val=50;
handles.xs=handles.z;
handles.ys=handles.y;
guidata(hObject,handles);
    cla;
if handles.segment==100
	img=label2rgb(squeeze(handles.Lrgb_fv(50,:,:)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(handles.I(50,:,:)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);
else
    axes(handles.axes1);
    imshow(squeeze(handles.final(1,:,:)));axis off;axis image;
end


% --------------------------------------------------------------------
function top_view_Callback(hObject, eventdata, handles)
% hObject    handle to top_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fv=0;
handles.tv=1;
handles.sv=0;
handles.slider_val=50;
handles.xs=handles.x;
handles.ys=handles.z;
guidata(hObject,handles);
    cla;
if handles.segment==100
	img=label2rgb(squeeze(handles.Lrgb_tv(:,50,:)),'jet','w','shuffle');
	axes(handles.axes1);
	imshow(squeeze(handles.I(:,50,:)));axis off;axis image;
	hold on;
	himage=imshow(img);axis off;axis image;
	set(himage,'Alphadata',0.3);
else
    axes(handles.axes1);
    imshow(squeeze(handles.final(:,1,:)));axis off;axis image;
end

% --------------------------------------------------------------------
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
segment=0;
handles.segment=double(segment);
guidata(hObject,handles);
p=round(handles.slider_val);
cla;
if handles.sv==1
    axes(handles.axes1);
    imshow(squeeze(handles.I(:,:,p)));
elseif handles.tv==1
    axes(handles.axes1);
    imshow(squeeze(handles.I(:,p,:)));
elseif handles.fv==1
    axes(handles.axes1);
    imshow(squeeze(handles.I(p,:,:)));
end


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function getting_started_Callback(hObject, eventdata, handles)
% hObject    handle to getting_started (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox('Press "load" to initiate process; use "clear" to remove watershed-mapping','Help-getting started');


function important_Callback(hObject, eventdata, handles)
% hObject    handle to important (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function navigation_Callback(hObject, eventdata, handles)
% hObject    handle to navigation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox('use "viewpoint" to view images in different planes & "scroll bar" to move along slices ','help-navigation');


% --------------------------------------------------------------------
function manual_segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to manual_segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox('the quality of segmentation depends on inital selection of area; NOTE-selecting region in one plane is sufficient to segment in 3D, different planes could be viewed using "viewpoint"','help-manual segmentation');

% --------------------------------------------------------------------
function automatic_segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to automatic_segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox('quality of segmentation depends on features of current image-slice displayed; NOTE-selecting region in one plane is sufficient to segment in 3D, different planes could be viewed using "viewpoint";','help-automatic segmentation');
