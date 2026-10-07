pausetime=10; % Pause during scan

%% Starting the laser scan and waiting for it to finish. Scan velocity MUST be set by hand
% Laser scan
las_chk=NP_USB_reperror_v2(NP_USB.Write(USBADDR, 'OUTP:SCAN:START'),'Scan Start');
if verboseflag
    display('*************************')
    display('Waiting for laser scan')
    display('*************************')
end
% waiting for the scan to finish
pause(pausetime);









