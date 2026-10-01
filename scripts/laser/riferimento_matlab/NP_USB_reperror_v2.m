function errc=NP_USB_reperror_v2(errorcode, optext)
% NP_USB_REPERROR report interpretation of Newport USB firmware error code
% Part of the Newport USB device Matlab code
% Adriaan Taal, Electrical Engineering - Columbia University

if errorcode  
    fprintf([optext ' operation failed \n'])
else
    fprintf([optext ' operation correctly executed \n'])
end

errc=errorcode;



