function varargout = arayuz(varargin)
% ARAYUZ MATLAB code for arayuz.fig
%      ARAYUZ, by itself, creates a new ARAYUZ or raises the existing
%      singleton*.
%
%      H = ARAYUZ returns the handle to a new ARAYUZ or the handle to
%      the existing singleton*.
%
%      ARAYUZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARAYUZ.M with the given input arguments.
%
%      ARAYUZ('Property','Value',...) creates a new ARAYUZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before arayuz_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to arayuz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help arayuz

% Last Modified by GUIDE v2.5 17-Jan-2017 01:06:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arayuz_OpeningFcn, ...
                   'gui_OutputFcn',  @arayuz_OutputFcn, ...
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


% --- Executes just before arayuz is made visible.
function arayuz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to arayuz (see VARARGIN)

% Choose default command line output for arayuz
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes arayuz wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = arayuz_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[y_original,Fs] = wavread('faruk');
canal_izquierdo_original=y_original(:,1);
t=linspace(0,length(canal_izquierdo_original)/Fs,length(canal_izquierdo_original));
plot(handles.axes1,t,canal_izquierdo_original);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[y_original,Fs,format] = wavread('faruk');
canal_izquierdo_original=y_original(:,1);

canal_izquierdo_original_reversa=flipud (canal_izquierdo_original);
t=linspace(0,length(canal_izquierdo_original)/Fs,length(canal_izquierdo_original));
plot(handles.axes2,t,canal_izquierdo_original_reversa);
%title(handles.axes2,'Fig 2. Ters Okuma')


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[y_original,Fs,format] = wavread('faruk');
size(y_original)   %2 canales (estereo): 569039 filas por 2columnas

canal_izquierdo_original=y_original(:,1); %
size(canal_izquierdo_original) %569039 filas por 1 columnas

sound(y_original,Fs)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[y_original,Fs,format] = wavread('faruk');
size(y_original)   %2 canales (estereo): 569039 filas por 2columnas

canal_izquierdo_original=y_original(:,1); %
size(canal_izquierdo_original) %569039 filas por 1 columnas

canal_izquierdo_original_reversa=flipud (canal_izquierdo_original);
sound(canal_izquierdo_original_reversa,Fs)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[speech,fs,nbits]=wavread('faruk.wav');  %Bölütlenecek ses dosyasý okunuyor.
window_ms = 100;  % Örnekleme periyodu belirleniyor
threshold = 1;     % Sesli bölge eþik deðeri belirleniyor
k_s=0;

% alinacak parcalarin sample sayisi
window = window_ms*fs/1000;  %alinacak parcalarin sample sayisi hesaplanýyor
speech2 = speech(1:(length(speech) - mod(length(speech),window)),1); % sinyal uzunlugu kontrolu
samples = reshape(speech2,window,length(speech2)/window); % yapýlýyor.

energies = sqrt(sum(samples.*samples))'; % Mühendislik öðrencilerinin aþina olduðu bir formül ;
%sinyalin samplelarýnýn karelerinin toplamý enerji deðerini veriyordu. Burada da bu iþlem yapýlarak
%energies dizisine sesin enerji deðerleri atýldý.
vuv = energies > threshold;% Bu satýrda eþik deðerimizi geçen örneklenmiþ aralýklarý 1 , eþik deðerin
%altýnda kalan satýrlarý 0 olarak atýyoruz.


%plot(energies);% komutu ile hesapladýðýmýz enerji dizisini çizdirip sesli  sessiz bölgeleri görme imkaný bulabiliyoruz. Kaydettiðim sesde 4 kelime bulunmakta. Bu grafikte bunu görebilmekteyiz.





%plot(vuv) ; %komutu ile eþik deðerlerini belirlediðimiz diziyi çizdirdik. Bu grafikte sesli ve sessiz bölgeler net bir biçimde tespit edilmiþ durumda. Þimdi basit bir algoritma geliþtirerek seslerimizi bölütleyip 4 yeni ses dosyasý halinde kaydedelim.



% Bölütlemede kullanýlan algoritma ; seslerin hangi örnekleme aralýðýnda baþlayýp hangi örnekleme aralýðýnda bittiðini belirleme üzerine çalýþmaktadýr.

z=length(energies);    %energies dizisinin uzunluðu belirleniyor.
a=1;
t=1;
%Bu döngüde  1 den sonra 0 gelen kaç örnekleme olduðunu tespit ediyoruz. Bu þekilde kaydýmýzda kaç kelime olduðunu anlýyoruz.Daha sonra bas ve son olmak üzere iki tane dizi tanýmlýyoruz.
        for j=1:z;
            if((vuv(j)==1) && (vuv(j+1)==0));
                k_s=k_s+1;
            end  
        end
bas=zeros(1,k_s);
son=zeros(1,k_s); 

% Bu döngüde 0 dan sonra bir gelen örneklerin(sample) yeri tespit ediliyor ve örnekleme aralýðý deðerimiz(window_ms) olan 4800 ile çarpýlarak kelimenin baþlangýç noktasý belirleniyor. Burada bas(i) i 'inci kelimenin baþlama noktasýný belirtiyor.
        for i=1:k_s
            for k=a:z
                if((vuv(k)==0) && (vuv(k+1)==1)); break;
                end
            end
            a=k+1;
            bas(i)=k*4400;
        end
%Bu döngüde 1 den sonra 0 gelen örneklerin (sample) yeri tespit ediliyor ve örnekleme aralýðý olan window_ms deðerimizle çarpýlarak kelimenin biriþ noktasý belirleniyor. Burada son(e)  e' inci kelimenin bitiþ noktasýný belirtiyor.

        for e=1:k_s
            for r=t:z
                if((vuv(r)==1) && (vuv(r+1)==0)); break;
                end
            end
            t=r+1;
            son(e)=r*4400;
        end

% Bölütlenen kelimeleri her birini yeni bir ses dosyasý halinde yazdýrýyoruz. Bu yazdýrma iþlemini detaylý bir þekilde anlatmayacaðým. Daha önceki yazýlarýmda buna detaylýca deðinmiþtim.
    k1=bas(1):son(1);
    plot(handles.axes3,speech(k1))
    %wavwrite(speech(k1),fs,'d1.wav');
    k2=bas(2):son(2);
    plot(handles.axes4,speech(k2))
    
    %wavwrite(speech(k2),fs,'d2.wav');
    k3=bas(3):son(3);
    plot(handles.axes5,speech(k3))
    %wavwrite(speech(k3),fs,'d3.wav');
    %k4=bas(4):son(4);
    %wavwrite(speech(k4),fs,'d4.wav');


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[speech,fs,nbits]=wavread('faruk.wav');  %Bölütlenecek ses dosyasý okunuyor.
window_ms = 100;  % Örnekleme periyodu belirleniyor
threshold = 1;     % Sesli bölge eþik deðeri belirleniyor
k_s=0;

% alinacak parcalarin sample sayisi
window = window_ms*fs/1000;  %alinacak parcalarin sample sayisi hesaplanýyor
speech2 = speech(1:(length(speech) - mod(length(speech),window)),1); % sinyal uzunlugu kontrolu
samples = reshape(speech2,window,length(speech2)/window); % yapýlýyor.

energies = sqrt(sum(samples.*samples))'; % Mühendislik öðrencilerinin aþina olduðu bir formül ;
%sinyalin samplelarýnýn karelerinin toplamý enerji deðerini veriyordu. Burada da bu iþlem yapýlarak
%energies dizisine sesin enerji deðerleri atýldý.
vuv = energies > threshold;% Bu satýrda eþik deðerimizi geçen örneklenmiþ aralýklarý 1 , eþik deðerin
%altýnda kalan satýrlarý 0 olarak atýyoruz.


%plot(energies);% komutu ile hesapladýðýmýz enerji dizisini çizdirip sesli  sessiz bölgeleri görme imkaný bulabiliyoruz. Kaydettiðim sesde 4 kelime bulunmakta. Bu grafikte bunu görebilmekteyiz.





%plot(vuv) ; %komutu ile eþik deðerlerini belirlediðimiz diziyi çizdirdik. Bu grafikte sesli ve sessiz bölgeler net bir biçimde tespit edilmiþ durumda. Þimdi basit bir algoritma geliþtirerek seslerimizi bölütleyip 4 yeni ses dosyasý halinde kaydedelim.



% Bölütlemede kullanýlan algoritma ; seslerin hangi örnekleme aralýðýnda baþlayýp hangi örnekleme aralýðýnda bittiðini belirleme üzerine çalýþmaktadýr.

z=length(energies);    %energies dizisinin uzunluðu belirleniyor.
a=1;
t=1;
%Bu döngüde  1 den sonra 0 gelen kaç örnekleme olduðunu tespit ediyoruz. Bu þekilde kaydýmýzda kaç kelime olduðunu anlýyoruz.Daha sonra bas ve son olmak üzere iki tane dizi tanýmlýyoruz.
        for j=1:z;
            if((vuv(j)==1) && (vuv(j+1)==0));
                k_s=k_s+1;
            end  
        end
bas=zeros(1,k_s);
son=zeros(1,k_s); 

% Bu döngüde 0 dan sonra bir gelen örneklerin(sample) yeri tespit ediliyor ve örnekleme aralýðý deðerimiz(window_ms) olan 4800 ile çarpýlarak kelimenin baþlangýç noktasý belirleniyor. Burada bas(i) i 'inci kelimenin baþlama noktasýný belirtiyor.
        for i=1:k_s
            for k=a:z
                if((vuv(k)==0) && (vuv(k+1)==1)); break;
                end
            end
            a=k+1;
            bas(i)=k*4400;
        end
%Bu döngüde 1 den sonra 0 gelen örneklerin (sample) yeri tespit ediliyor ve örnekleme aralýðý olan window_ms deðerimizle çarpýlarak kelimenin biriþ noktasý belirleniyor. Burada son(e)  e' inci kelimenin bitiþ noktasýný belirtiyor.

        for e=1:k_s
            for r=t:z
                if((vuv(r)==1) && (vuv(r+1)==0)); break;
                end
            end
            t=r+1;
            son(e)=r*4400;
        end

% Bölütlenen kelimeleri her birini yeni bir ses dosyasý halinde yazdýrýyoruz. Bu yazdýrma iþlemini detaylý bir þekilde anlatmayacaðým. Daha önceki yazýlarýmda buna detaylýca deðinmiþtim.
    k1=bas(1):son(1);
    sound(speech(k1),fs);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[speech,fs,nbits]=wavread('faruk.wav');  %Bölütlenecek ses dosyasý okunuyor.
window_ms = 100;  % Örnekleme periyodu belirleniyor
threshold = 1;     % Sesli bölge eþik deðeri belirleniyor
k_s=0;

% alinacak parcalarin sample sayisi
window = window_ms*fs/1000;  %alinacak parcalarin sample sayisi hesaplanýyor
speech2 = speech(1:(length(speech) - mod(length(speech),window)),1); % sinyal uzunlugu kontrolu
samples = reshape(speech2,window,length(speech2)/window); % yapýlýyor.

energies = sqrt(sum(samples.*samples))'; % Mühendislik öðrencilerinin aþina olduðu bir formül ;
%sinyalin samplelarýnýn karelerinin toplamý enerji deðerini veriyordu. Burada da bu iþlem yapýlarak
%energies dizisine sesin enerji deðerleri atýldý.
vuv = energies > threshold;% Bu satýrda eþik deðerimizi geçen örneklenmiþ aralýklarý 1 , eþik deðerin
%altýnda kalan satýrlarý 0 olarak atýyoruz.


%plot(energies);% komutu ile hesapladýðýmýz enerji dizisini çizdirip sesli  sessiz bölgeleri görme imkaný bulabiliyoruz. Kaydettiðim sesde 4 kelime bulunmakta. Bu grafikte bunu görebilmekteyiz.





%plot(vuv) ; %komutu ile eþik deðerlerini belirlediðimiz diziyi çizdirdik. Bu grafikte sesli ve sessiz bölgeler net bir biçimde tespit edilmiþ durumda. Þimdi basit bir algoritma geliþtirerek seslerimizi bölütleyip 4 yeni ses dosyasý halinde kaydedelim.



% Bölütlemede kullanýlan algoritma ; seslerin hangi örnekleme aralýðýnda baþlayýp hangi örnekleme aralýðýnda bittiðini belirleme üzerine çalýþmaktadýr.

z=length(energies);    %energies dizisinin uzunluðu belirleniyor.
a=1;
t=1;
%Bu döngüde  1 den sonra 0 gelen kaç örnekleme olduðunu tespit ediyoruz. Bu þekilde kaydýmýzda kaç kelime olduðunu anlýyoruz.Daha sonra bas ve son olmak üzere iki tane dizi tanýmlýyoruz.
        for j=1:z;
            if((vuv(j)==1) && (vuv(j+1)==0));
                k_s=k_s+1;
            end  
        end
bas=zeros(1,k_s);
son=zeros(1,k_s); 

% Bu döngüde 0 dan sonra bir gelen örneklerin(sample) yeri tespit ediliyor ve örnekleme aralýðý deðerimiz(window_ms) olan 4800 ile çarpýlarak kelimenin baþlangýç noktasý belirleniyor. Burada bas(i) i 'inci kelimenin baþlama noktasýný belirtiyor.
        for i=1:k_s
            for k=a:z
                if((vuv(k)==0) && (vuv(k+1)==1)); break;
                end
            end
            a=k+1;
            bas(i)=k*4400;
        end
%Bu döngüde 1 den sonra 0 gelen örneklerin (sample) yeri tespit ediliyor ve örnekleme aralýðý olan window_ms deðerimizle çarpýlarak kelimenin biriþ noktasý belirleniyor. Burada son(e)  e' inci kelimenin bitiþ noktasýný belirtiyor.

        for e=1:k_s
            for r=t:z
                if((vuv(r)==1) && (vuv(r+1)==0)); break;
                end
            end
            t=r+1;
            son(e)=r*4400;
        end

% Bölütlenen kelimeleri her birini yeni bir ses dosyasý halinde yazdýrýyoruz. Bu yazdýrma iþlemini detaylý bir þekilde anlatmayacaðým. Daha önceki yazýlarýmda buna detaylýca deðinmiþtim.
    
    k2=bas(2):son(2);
    sound(speech(k2),fs);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[speech,fs,nbits]=wavread('faruk.wav');  %Bölütlenecek ses dosyasý okunuyor.
window_ms = 100;  % Örnekleme periyodu belirleniyor
threshold = 1;     % Sesli bölge eþik deðeri belirleniyor
k_s=0;

% alinacak parcalarin sample sayisi
window = window_ms*fs/1000;  %alinacak parcalarin sample sayisi hesaplanýyor
speech2 = speech(1:(length(speech) - mod(length(speech),window)),1); % sinyal uzunlugu kontrolu
samples = reshape(speech2,window,length(speech2)/window); % yapýlýyor.

energies = sqrt(sum(samples.*samples))'; % Mühendislik öðrencilerinin aþina olduðu bir formül ;
%sinyalin samplelarýnýn karelerinin toplamý enerji deðerini veriyordu. Burada da bu iþlem yapýlarak
%energies dizisine sesin enerji deðerleri atýldý.
vuv = energies > threshold;% Bu satýrda eþik deðerimizi geçen örneklenmiþ aralýklarý 1 , eþik deðerin
%altýnda kalan satýrlarý 0 olarak atýyoruz.


%plot(energies);% komutu ile hesapladýðýmýz enerji dizisini çizdirip sesli  sessiz bölgeleri görme imkaný bulabiliyoruz. Kaydettiðim sesde 4 kelime bulunmakta. Bu grafikte bunu görebilmekteyiz.





%plot(vuv) ; %komutu ile eþik deðerlerini belirlediðimiz diziyi çizdirdik. Bu grafikte sesli ve sessiz bölgeler net bir biçimde tespit edilmiþ durumda. Þimdi basit bir algoritma geliþtirerek seslerimizi bölütleyip 4 yeni ses dosyasý halinde kaydedelim.



% Bölütlemede kullanýlan algoritma ; seslerin hangi örnekleme aralýðýnda baþlayýp hangi örnekleme aralýðýnda bittiðini belirleme üzerine çalýþmaktadýr.

z=length(energies);    %energies dizisinin uzunluðu belirleniyor.
a=1;
t=1;
%Bu döngüde  1 den sonra 0 gelen kaç örnekleme olduðunu tespit ediyoruz. Bu þekilde kaydýmýzda kaç kelime olduðunu anlýyoruz.Daha sonra bas ve son olmak üzere iki tane dizi tanýmlýyoruz.
        for j=1:z;
            if((vuv(j)==1) && (vuv(j+1)==0));
                k_s=k_s+1;
            end  
        end
bas=zeros(1,k_s);
son=zeros(1,k_s); 

% Bu döngüde 0 dan sonra bir gelen örneklerin(sample) yeri tespit ediliyor ve örnekleme aralýðý deðerimiz(window_ms) olan 4800 ile çarpýlarak kelimenin baþlangýç noktasý belirleniyor. Burada bas(i) i 'inci kelimenin baþlama noktasýný belirtiyor.
        for i=1:k_s
            for k=a:z
                if((vuv(k)==0) && (vuv(k+1)==1)); break;
                end
            end
            a=k+1;
            bas(i)=k*4400;
        end
%Bu döngüde 1 den sonra 0 gelen örneklerin (sample) yeri tespit ediliyor ve örnekleme aralýðý olan window_ms deðerimizle çarpýlarak kelimenin biriþ noktasý belirleniyor. Burada son(e)  e' inci kelimenin bitiþ noktasýný belirtiyor.

        for e=1:k_s
            for r=t:z
                if((vuv(r)==1) && (vuv(r+1)==0)); break;
                end
            end
            t=r+1;
            son(e)=r*4400;
        end

% Bölütlenen kelimeleri her birini yeni bir ses dosyasý halinde yazdýrýyoruz. Bu yazdýrma iþlemini detaylý bir þekilde anlatmayacaðým. Daha önceki yazýlarýmda buna detaylýca deðinmiþtim.
    
    k3=bas(3):son(3);
    sound(speech(k3),fs);


% --- Executes during object creation, after setting all properties.
function uipanel5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ses,Fs] = wavread('faruk');
duzSes=ses(30800:33800,1);
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes6,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(34700:37800,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes22,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(39900:43000,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes25,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ses,Fs] = wavread('faruk');
duzSes=ses(42500:45300,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes26,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');
% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ses,Fs] = wavread('faruk');
duzSes=ses(49100:54000,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes27,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(66000:67100,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes28,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk'); 
duzSes=ses(67800:71000,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes29,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(74500:76500,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes30,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');
% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(76500:80000,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes31,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(80700:85000,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes32,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');
% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ses,Fs] = wavread('faruk');
duzSes=ses(89100:92600,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes33,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');
% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(92350:94800,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes34,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');
% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(96700:99800,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes35,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');

% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(99500:109300,1);
sound(duzSes, Fs)
t=linspace(0,length(duzSes)/Fs,length(duzSes));
plot(handles.axes36,t,duzSes)

tersSes=flipud(duzSes);
sound(duzSes, Fs)

set(handles.figure1,'HandleVisibility','off');
    close all
set(handles.figure1,'HandleVisibility','on');

ydftd = fft(duzSes);
ydftt = fft(tersSes);
% I'll assume y has even length
ydftd = ydftd(1:length(duzSes)/2+1);
ydftt = ydftt(1:length(tersSes)/2+1);
% create a frequency vector
freqd = 0:Fs/length(duzSes):Fs/2;
freqt = 0:Fs/length(tersSes):Fs/2;
vd=duzSes(:,1);
vt=tersSes(:,1);
figure;
subplot(3,2,1);
plot(vd);
title('Signal')
xlabel('Time');

subplot(3,2,2);
plot(vt);
title('Reverse of Signal')
xlabel('Time');
% plot magnitude
subplot(3,2,3);
plot(freqd,abs(ydftd));
title('Magnitude of FFT')
xlabel('Hz');

subplot(3,2,4);
plot(freqt,abs(ydftt));
title('Magnitude of Reverse FFT')
xlabel('Hz');
% plot phase
subplot(3,2,5);
plot(freqd,unwrap(angle(ydftd)));
title('Phase of FFT')
xlabel('Hz');

subplot(3,2,6);
plot(freqt,unwrap(angle(ydftt)));
title('Phase of Reverse FFT')
xlabel('Hz');


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
J = imrotate(matlabImage,-45,'bilinear');
image(J)
axis off
axis image


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
J = imrotate(matlabImage,45,'bilinear');
image(J)
axis off
axis image


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
J = imrotate(matlabImage,-90,'bilinear');
image(J)
axis off
axis image


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
J = imrotate(matlabImage,90,'bilinear');
image(J)
axis off
axis image

% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
J = imrotate(matlabImage,180,'bilinear');
image(J)
axis off
axis image


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes39)
matlabImage = imread('vesika.jpg');
image(matlabImage)
axis off
axis image

% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
imgMirror = flipdim(matlabImage,2);
%J = imrotate(matlabImage,180,'bilinear');
image(imgMirror)
axis off
axis image


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
imgMirror = flipdim(matlabImage,1);
%J = imrotate(matlabImage,180,'bilinear');
image(imgMirror)
axis off
axis image


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes40)
matlabImage = imread('vesika.jpg');
J = rgb2gray(matlabImage);
imshow(J)
axis off
axis image


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[y,fs] = wavread('faruk.wav');
v=y(:,1);
ydft = fft(v);
ydft = ydft(1:length(v)/2+1);
freq = 0:fs/length(v):fs/2;


vTers=flipud(v);

Tersfft = fft(vTers);
Tersfft = Tersfft(1:length(vTers)/2+1);
freqTers = 0:fs/length(vTers):fs/2;



figure;

subplot(3,2,1);
plot(v);
title('Signal')
xlabel('Zaman');

subplot(3,2,2);
plot(vTers);
title('Reverse of Signal')
xlabel('Zaman');

subplot(3,2,3);
plot(freq,abs(ydft));
title('FFT')
xlabel('Frequency(Hz)');

subplot(3,2,4);
plot(freqTers,abs(Tersfft));
title('FFT of Reverse Signal')
xlabel('Frequency(Hz)');

subplot(3,2,5);
plot(freq,unwrap(angle(ydft)));
title('Phase FFT')
xlabel('Frequency(Hz)');

subplot(3,2,6);
plot(freqTers,unwrap(angle(Tersfft)));
title('Phase FFT of Reverse Signal')
xlabel('Frequency(Hz)');


% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[y,fs] = wavread('faruk.wav');
v=y(:,1);
asd=length(y)/fs;
dd=length(y);
t=linspace(0,length(y)/fs,length(y));
ters=v';

fftStd = fft(v);
fftStd = fftStd(1:length(v)/2+1);
freq = 0:fs/length(v):fs/2;

figure;

subplot(4,4,2);
plot(t,v);
title('Signal')

subplot(4,4,3);
plot(freq,abs(fftStd));
title('Sinyal FFT')

subplot(4,4,4);
plot(freq,unwrap(angle(fftStd)));
title('Sinyal FFT Faz')

f=5;
kos=cos(2*pi*t*f);
subplot(4,4,5)
plot(t,kos)
title('cos (f=5)')

yeni=ters.*cos(2*pi*t*f);
subplot(4,4,6)
plot(t,yeni)
title('Cos x Signal')

fftStd = fft(yeni);
fftStd = fftStd(1:length(yeni)/2+1);

subplot(4,4,7)
plot(freq,abs(fftStd))
title('Cos x Signal FFT')

subplot(4,4,8)
plot(freq,unwrap(angle(fftStd)))
title('Cos x Signal FFT')

f=50;
kos=cos(2*pi*t*f);
subplot(4,4,9)
plot(t,kos)
title('cos (f=50)')


yeni=ters.*cos(2*pi*t*f);
subplot(4,4,10)
plot(t,yeni)
title('Cos x Signal')

fftStd = fft(yeni);
fftStd = fftStd(1:length(yeni)/2+1);

subplot(4,4,11)
plot(freq,abs(fftStd))
title('Cos x Signal FFT')

subplot(4,4,12)
plot(freq,unwrap(angle(fftStd)))
title('Cos x Signal FFT')

f=5000;
kos=cos(2*pi*t*f);
subplot(4,4,13)
plot(t,kos)
title('cos (f=5000)')

yeni=ters.*cos(2*pi*t*f);
subplot(4,4,14)
plot(t,yeni)
title('Cos x Signal')

fftStd = fft(yeni);
fftStd = fftStd(1:length(yeni)/2+1);

subplot(4,4,15)
plot(freq,abs(fftStd))
title('Cos x Signal FFT')

subplot(4,4,16)
plot(freq,unwrap(angle(fftStd)))
title('Cos x Signal FFT')


% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[normalFrekans,fs] = wavread('faruk.wav');
nF=normalFrekans(:,1);
nt=[1/fs:1/fs:length(normalFrekans)/fs];

wavwrite(nF,fs/2,'dusuk');
wavwrite(nF,2*fs,'yuksek');

[dusukFrekans,dfs] = wavread('dusuk.wav');
dF=dusukFrekans(:,1);
dt=[1/(fs/2):1/(fs/2):length(dusukFrekans)/(fs/2)];

dfft = fft(dF);
dfft = dfft(1:length(dF)/2+1);
dfreq = 0:dfs/length(dF):dfs/2;

[yuksekFrekans,yfs] = wavread('yuksek.wav');
yF=yuksekFrekans(:,1);
yt=[1/(2*fs):1/(2*fs):length(yuksekFrekans)/(2*fs)];

yfft = fft(yF);
yfft = yfft(1:length(yF)/2+1);
yfreq = 0:yfs/length(yF):yfs/2;

figure;
subplot(4,2,[1 2]);
plot(nt,nF);
title('Sinyal')

subplot(4,2,3);
plot(dt,dF);
title('Düþük Frekans')

subplot(4,2,4);
plot(yt,yF);
title('Yüksek Frekans')

subplot(4,2,5);
plot(dfreq,abs(dfft));
title('Düþük Frekans FFT')

subplot(4,2,6);
plot(yfreq,abs(yfft));
title('Yüksek Frekans FFT')

subplot(4,2,7);
plot(dfreq,unwrap(angle(dfft)));
title('Düþük Frekans FFT Faz')
subplot(4,2,8);
plot(yfreq,unwrap(angle(yfft)));
title('Yüksek Frekans FFT Faz')


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[yuksekFrekans,yfs] = wavread('yuksek.wav');
yF=yuksekFrekans(:,1);
sound(yF,yfs);

% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dusukFrekans,dfs] = wavread('dusuk.wav');
dF=dusukFrekans(:,1);
sound(dF,dfs);


% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ses,Fs] = wavread('faruk');
duzSes=ses(30800:33800,1);
t=linspace(0,length(duzSes)/Fs,length(duzSes));

Nfft = fft(duzSes);
Nfft = Nfft(1:length(duzSes)/2+1);
Nfreq = 0:Fs/length(duzSes):Fs/2;
% tersSes=flipud(duzSes);

sonSes=resample(duzSes,2,1);
st=linspace(0,length(sonSes)/Fs,length(sonSes));
% sound(sonSes,Fs);

Sfft = fft(sonSes);
Sfft = Sfft(1:length(sonSes)/2+1);
Sfreq = 0:Fs/length(sonSes):Fs/2;

son2Ses=resample(duzSes,1,2);
s2t=linspace(0,length(son2Ses)/Fs,length(son2Ses));
% sound(son2Ses,Fs);

S2fft = fft(son2Ses);
S2fft = S2fft(1:length(son2Ses)/2+1);
S2freq = 0:Fs/length(son2Ses):Fs/2;

figure;
subplot(3,3,1);
plot(t,duzSes)
title('Normal A Signal');

subplot(3,3,2);
plot(st,sonSes)
title('More Sample for A Signal');

subplot(3,3,3);
plot(s2t,son2Ses)
title('Less Sample for A Signal');

subplot(3,3,4);
plot(Nfreq,abs(Nfft));
title('Normal A Signal FFT');

subplot(3,3,5);
plot(Sfreq,abs(Sfft));
title('More Sample for A Signal FFT');

subplot(3,3,6);
plot(S2freq,abs(S2fft));
title('More Sample for A Signal FFT');

subplot(3,3,7);
plot(Nfreq,unwrap(angle(Nfft)));
title('Normal A Signal FFT Phase');

subplot(3,3,8);
plot(Sfreq,unwrap(angle(Sfft)));
title('Normal A Signal FFT Phase');

subplot(3,3,9);
plot(S2freq,unwrap(angle(S2fft)));
title('Normal A Signal FFT Phase');





% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
