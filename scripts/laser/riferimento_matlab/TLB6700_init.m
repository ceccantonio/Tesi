%% Script for initializing the TLB6700 laser. 
% In this version, reperror has been modified to ignore the USBDUPLICATEADDRESS error


deviceID = hex2dec('100A');
dll_path=[pwd,'\TLB6700'];

addpath(dll_path);

NPasm = NET.addAssembly([dll_path '\UsbDllWrap.dll']);


%Get a handle on the USB class
NPASMtype = NPasm.AssemblyHandle.GetType('Newport.USBComm.USB');
%launch the class USB, it constructs and allows to use functions in USB.h
NP_USB = System.Activator.CreateInstance(NPASMtype);

%Open the USB device
NP_USB_reperror(NP_USB.OpenDevices(deviceID),'DeviceOpen');

%The Query method sends the passed in command string to the specified device and reads the response data.
querydata = System.Text.StringBuilder(64);
las_chk=NP_USB_reperror_v2(NP_USB.Query(USBADDR, '*IDN?', querydata),'Query');
devInfo = char(ToString(querydata));
fprintf(['Device attached is ' devInfo '\n']);

% Check the instrument output status (ON/OFF)
NP_USB_reperror(NP_USB.Query(USBADDR, 'OUTP:STAT?', querydata),'Query');
devState = char(ToString(querydata));
fprintf(['Device state is ' devState '\n']);
