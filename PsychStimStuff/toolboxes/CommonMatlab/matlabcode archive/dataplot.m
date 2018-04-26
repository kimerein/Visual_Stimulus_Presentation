%Display Data
figure(10)
subplot(2,1,1)
plot(intraCell_data(:,20,1))
subplot(2,1,2)
plot(extraCell_data(:,20,1))

figure(10)
subplot(2,1,1)
plot(intraCell_data(:,20,2))
subplot(2,1,2)
plot(extraCell_data(:,20,2))


Scalar =  ones(size(extracted_spikes,1),1)*extracted_spikes(63,:,4);
figure(5)
subplot(2,1,1)
plot(extracted_spikes(:,:,1))
subplot(2,1,2)
plot(extracted_spikes(:,:,3)./Scalar )
plot(extracted_spikes(:,:,3) )



figure(5)
subplot(2,1,1)
plot(extracted_spikes(:,find(extracted_spikes(4,:,2) == 1,1)))
plot(extracted_spikes(:,find(extracted_spikes(4,:,2) == 2,1)))
subplot(2,1,2)
plot(extracted_spikes(:,find(extracted_spikes(4,:,4) == 1,3)))
plot(extracted_spikes(:,find(extracted_spikes(4,:,4) == 2,3)))




