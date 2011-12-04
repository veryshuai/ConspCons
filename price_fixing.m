function fix = price_fixing(p)
%this function linearly extends prices (according to 4yr trend) where I don't have data, and
%normalizes to 2010  

ind = zeros(1,size(p,2));
fix = p;
for k = 1:size(p,2);
ind(k) = find(p(:,k),1);
if ind(k) ~= 1
   b = polyfit((1:4)',p(ind(k):ind(k)+3,k),1); 
   for m = 1:ind(k)-1
      fix(m,k) = p(ind(k),k)-(ind(k)-m)*b(1); 
   end
end
fix(:,k) = fix(:,k)/fix(end,k);
end

end