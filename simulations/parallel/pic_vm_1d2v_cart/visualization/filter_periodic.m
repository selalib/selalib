% Periodic binomial filter
% in: input data
% out2: filtered data
%
% author: Katharina Kormann, 2019
function out2 = filter_periodic(in)

n = size(in);
out = zeros(n(1),n(2));

out(1,:) = (in(n(1),:)+in(2,:))*0.25+in(1,:)*0.5;
for i=2:n(1)-1
    out(i,:) = (in(i-1,:)+in(i+1,:))*0.25+in(i,:)*0.5;
end
out(n(1),:) = (in(n(1)-1,:)+in(1,:))*0.25+in(n(1),:)*0.5;

out2 = zeros(n(1),n(2));
 
out2(:,1) = (out(:,n(2))+out(:,2))*0.25+out(:,1)*0.5;
for i=2:n(2)-1
    out2(:,i) = (out(:,i-1)+out(:,i+1))*0.25+out(:,i)*0.5;
end
out2(:,n(2)) = (out(:,n(2)-1)+out(:,1))*0.25+out(:,n(2))*0.5;
