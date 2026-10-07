function vai = NP_USB_reperror(errorcode, optext)
% NP_USB_REPERROR report interpretation of Newport USB firmware error code
% Part of the Newport USB device Matlab code
% Adriaan Taal, Electrical Engineering - Columbia University

switch errorcode
    case 0
        fprintf([optext ' operation correctly executed \n'])
    case 1
        fprintf([optext ' operation correctly executed \n'])
        %fprintf([optext ' error USBDUPLICATEADDRESS, More than one device on the bus has the same device ID \n'])
    case -2
        fprintf([optext ' error USBADDRESSNOTFOUND, The device ID cannot be found among the open devices on the bus \n'])
    case -3
        fprintf([optext ' error USBINVALIDADDRESS, The device ID is outside the valid range of 0 - 31 \n'])
    otherwise
        fprintf([optext ' error Unknown... \n'])
end

vai = errorcode;

end

