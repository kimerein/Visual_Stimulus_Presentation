aa = unique(a(:,1));
for i =4:4 %length(aa)
    ind = find(a(:,1)==aa(i));
    b(round(a(ind,2)*1)) = a(ind,3); %% make into time array
end
ttt = xcov(b,2000,'coeff');
plot(ttt)
% plot(([1:size(ttt,1)]-((size(ttt,1)-1)/2 +1)),ttt)
% plot(b)