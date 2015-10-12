


function varargout = projectGUI(varargin)
% PROJECTGUI MATLAB code for projectGUI.fig
%      PROJECTGUI, by itself, creates a new PROJECTGUI or raises the existing
%      singleton*.
%
%      H = PROJECTGUI returns the handle to a new PROJECTGUI or the handle to
%      the existing singleton*.
%
%      PROJECTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECTGUI.M with the given input arguments.
%
%      PROJECTGUI('Property','Value',...) creates a new PROJECTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before projectGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to projectGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help projectGUI

% Last Modified by GUIDE v2.5 07-Jun-2013 13:48:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @projectGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @projectGUI_OutputFcn, ...
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
chkq=0;
handles.chk=chkq;% End initialization code - DO NOT EDIT
handles.chk1=chkq;

% --- Executes just before projectGUI is made visible.
function projectGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to projectGUI (see VARARGIN)

% Choose default command line output for projectGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projectGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = projectGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% load=uigetfile({'*.jpg;*.tif;*.bmp;*.png;*.gif','All Image Files';...
%           '*.*','All Files' },'Browse for Image',...
%           'D:\proimages\Edh0.jpg');
% handles.load=load;
[name surname]=uigetfile('*.*');
img=imread([surname name]);
% imtool(img);
imshow(img);
handles.obj=img;

guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[side1 maxArea1 minorr type img1]=fullfunc(handles.obj);
% fprintf(fid,'The study of non contrast CT scan of %d reveals that  %d type of bleeding ',handles.age,type);
imshow(img1(:,:,1));
winopen('C:\Program Files\MATLAB\R2012b\bin\ttttt.txt');
 handles.fidd = fopen('C:\Program Files\MATLAB\R2012b\bin\ttttt.txt', 'wt');
 fprintf(handles.fidd,'\n\n\n NAME----  %s  ',handles.name);
 fprintf(handles.fidd,'\n AGE----%d ',handles.age);
 if(handles.chk>=1)
     inj='HYPERTENSIVE';
 else
     inj='NORMOTENSIVE';
 end
 fprintf(handles.fidd,'\n PREVIOUS HEAD INJURY----%s ',inj);
 if(handles.chk1>=1)
     hy='YES';
 else
     hy='NO';
 end
 
 fprintf(handles.fidd,'\n HISTORY----%s ',hy);
 fprintf(handles.fidd,'\n GENDER----%s ',handles.gen);
  fprintf(handles.fidd,'\n\n\n\nThe non contrast study of CT scan of patient reavels  %s  bleed of thickness %d pixels and area %d pixels on the %s side of brain.',type,round(minorr),maxArea1,side1); 
 fprintf(handles.fidd,'\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\t\t\t\t\t\t\t\t\t\t\tDoctors Signature');
 
 
 
guidata(hObject, handles);





function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.age=   str2double(get(hObject,'String')) ;
age=handles.age
guidata(hObject, handles);
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
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
handles.chk= get(hObject,'Value')

guidata(hObject, handles);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
handles.chk1= get(hObject,'Value');

guidata(hObject, handles);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.name=get(hObject,'String')
 name1=handles.name;
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
handles.item= cellstr(get(hObject,'String')) 
% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 handles.cont1 = get(hObject,'Value');
 if(handles.cont1>=2)
     handles.gen='female';
 else
     handles.gen='male';%        contents{get(hObject,'Value')} returns selected item from listbox3
 end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function [side1 maxArea1 minorr type img1]=fullfunc(im)

% disp('here');

% side=z(1);
% bleed =z(2);
% maxArea1=z(3);
% minorr=z(4) 
% type=z(5);
% text=num2str(side);
% text=strcat(text,num2str(bleed));
% text=strcat(text,num2str(maxArea1));
% text=strcat(text,num2str(minorr));
% text=strcat(text,num2str(type));
% newimg=im>150;


im1=im;
image=im;
% image=imread('C:\Users\vishal\Desktop\test\18.jpg');
[m n] = size(image);
im2=rgb2gray(image);
q=1;

im3=double(im2);
im4=imfill(im2,'holes');
[m,n]=size(im2);
for i=1:m
    for j=1:n
        if( im4(i,j)<=230)
            z1(i,j)=0;
        else
            z1(i,j)=im4(i,j);
        end
        
    end
end
%  figure,imshow(im2);


for i=1:m
    for j=1:n
        im5(i,j)=(z1(i,j)-im2(i,j));
    end
end
%figure,imshow(im5);

for i=1:m
    for j=1:n
        im6(i,j)=255-im5(i,j);
    end
end

%figure,imshow(im6);
for i=1:m
    for j=1:n
        if( im6(i,j)>240)
            z(i,j)=0;
        else
            z(i,j)=im6(i,j);
        end
        
    end
end

z = medfilt2(z,[7 7]);
y=z;
%figure,imshow(uint8(y));

%figure,imshow(I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HISTOGRAM THRESHOLDING%%%%%%%%%%%%%%%%%

histthres = histogramthresfunc(im2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:m
    for j=1:n
        if(y(i,j)>=histthres && y(i,j)<=210)
            z2(i,j)=255;
        else
            z2(i,j)=0;
        end
        
    end
end

z2 = bwmorph(z2,'clean');


level =0.1;
bw = im2bw(z2,level);
%figure, imshow(bw);                       %3 steps to auto thres grey to binary image
cc = bwconncomp(bw);
numPixels = cellfun(@numel,cc.PixelIdxList);
stats = regionprops(cc, 'All');
[biggest,idx] = max(numPixels);
angle1=stats(idx).Orientation;
se = strel('line',3,angle1);
I3=bw;
for i=1:3
I2 = imdilate(I3,se);
I3=I2;
end
cc = bwconncomp(I3);
numPixels = cellfun(@numel,cc.PixelIdxList);
[biggest,idx] = max(numPixels);
BW2 = ismember(labelmatrix(cc), idx);
% figure,imshow(BW2);
title('bleed');
 img1=BW2;
%  img1=gray2ind(img1);
%  img1=ind2rgb(img1,[1 1 1]);
% for i=1:m
%     for j=1:n
%         if( BW2(i,j)>=1)
%                 img1(i,j)=[1 0 0];
% %         else
%             img1(i,j)=im1(i,j);
%         end
%         
%     end
% end
img1=double(img1);

CC = bwconncomp(BW2);
S = regionprops(CC,'all');
% figure,imshow(BW2);
title('bleed');

maxArea1 = max([S.Area]);
maxArea1=maxArea1-1;

if(maxArea1>200)
  load svm1.mat;
idx = find([S.Area] >= maxArea1);
side = sidefunc(BW2,S,idx);
if(side>=1 && side<2)
    side1='left';
elseif(side<=0)
    side1='right';
 end
minorr=S(idx).MinorAxisLength;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 contactarea = contactareafunc (im2,BW2,angle1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(contactarea>20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[solid,axisarea,Con,Ene,Homo,Cor,convexx]  = featurefunc(S,maxArea1,BW2,im2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pqr=[Con Homo Ene Cor solid axisarea];
Group = svmclassify(svmStruct,pqr,'Showplot',true);
if(Group>=1)
%     disp('SDH'); 
    type='Subdural Hemmorrhage';
    
    
else
%     disp('EDH');
    type='Epidural Hemorrhage';
end

 else
%      display('hypertensive bleed');
     type='Hypertensive';
 end
 else
%      display('no bleeding');
     type='No ';
     side=2
     minorr=0
     maxArea1=0
     img1=im;
     side1='either side';

end
 


% function edit5_Callback(hObject, eventdata, handles)
% gender=get(hObject,'String')
% guidata(hObject,handles);
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit5_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
