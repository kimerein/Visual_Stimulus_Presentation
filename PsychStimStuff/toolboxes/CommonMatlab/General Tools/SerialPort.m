s = serial('COM4','BaudRate',19200,'Parity','none','DataBits',8,'StopBit',1,'Terminator','CR','FlowControl','hardware')
fopen(s)
fprintf(s,'1AC?\n')
idn = fscanf(s)
fclose(s)
delete(s)
clear s
