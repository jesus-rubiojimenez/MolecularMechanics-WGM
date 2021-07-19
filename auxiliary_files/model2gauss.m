function [proplikelihood] = model2gauss(time,lshift,A1,A2,W1,W2,T1,T2,sigma0)
%model2gauss(time,A1,A2,W1,W2,T1,T2,sigma0)

temp = 0;
for i=1:length(time)
    temp=temp+(lshift(i)-A1.*exp(-(time(i)-T1).^2/(2.*W1.^2))+A2.*exp(-(time(i)-T2).^2/(2.*W2.^2))).^2;
end
proplikelihood = exp(-temp/(2*sigma0^2));

end

