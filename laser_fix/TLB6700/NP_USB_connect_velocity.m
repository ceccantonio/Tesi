clear all;
clc
USBADDR = 1;
deviceID = hex2dec('100A');
dll_path=pwd;

NPasm = NET.addAssembly([dll_path '\UsbDllWrap.dll']);

%Get a handle on the USB class
NPASMtype = NPasm.AssemblyHandle.GetType('Newport.USBComm.USB');
%launch the class USB, it constructs and allows to use functions in USB.h
NP_USB = System.Activator.CreateInstance(NPASMtype);

% %Open the USB device
% NP_USB_reperror(NP_USB.OpenDevices(deviceID),'DeviceOpen');
% 
% %Initialize Event handling
% NP_USB_reperror(NP_USB.EventInit(deviceID),'EventInit');

%Open the USB device
NP_USB_reperror(NP_USB.OpenDevices(deviceID),'DeviceOpen');

%The Query method sends the passed in command string to the specified device and reads the response data.
querydata = System.Text.StringBuilder(64);
NP_USB_reperror(NP_USB.Query(USBADDR, '*IDN?', querydata),'Query');
devInfo = char(ToString(querydata));
fprintf(['Device attached is ' devInfo '\n']);

%%

NP_USB_reperror(NP_USB.Query(USBADDR, 'OUTP:STAT?', querydata),'Query');
devState = char(ToString(querydata));
fprintf(['Device state is ' devState '\n']);

%%
%NP_USB_reperror(NP_USB.Write(USBADDR, 'OUTP:SCAN:START'));
%NP_USB.Write(USBADDR, 'OUTP:SCAN:START');

%%
% l_set=1535.0;
% 
% % NP_USB_reperror(NP_USB.Query(USBADDR, 'SENS:WAVE?', querydata),'Query');
% NP_USB.Query(USBADDR, 'SENS:WAVE?', querydata);
% lambda = str2double(char(ToString(querydata)));
% disp(['Initial Lambda = ', num2str(lambda)]);
% 
% NP_USB_reperror(NP_USB.Write(USBADDR, ['SOUR:WAVE ' num2str(l_set)]),['Wavelength set to ' num2str(l_set)]);
% tr=1;
% NP_USB_reperror(NP_USB.Write(USBADDR, ['OUTP:TRAC ' num2str(tr)]),['Track mode is ' num2str(tr)]);
% 
% disp('***********************************')
% % velocity 1 nm/s 
% % 1 s settling time
% 
% disp('waiting for lambda tracking')
% dl=abs(lambda-l_set);
% pause(dl+1);
% 
% disp('***********************************')
% 
% tr=0;
% NP_USB_reperror(NP_USB.Write(USBADDR, ['OUTP:TRAC ' num2str(tr)]),['Track mode is ' num2str(tr)]);
% 
% NP_USB_reperror(NP_USB.Query(USBADDR, 'SENS:WAVE?', querydata),'Query');
% Lambda = str2double(char(ToString(querydata)));

%% Start scan

NP_USB_reperror(NP_USB.Write(USBADDR, 'OUTP:SCAN:START'));

