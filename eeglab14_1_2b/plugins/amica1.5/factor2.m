function out = factor2(n)

f = factor(n);
if length(f) == 1
    f = [1 f];
end
numfact = length(f);

% first way
i = 1;
j = numfact;
r = [1 1];
k = 1;
while i < j
   if mod(k,2)
      r(1) = r(1)*f(i);
      r(2) = r(2)*f(j);
   else
      r(1) = r(1)*f(j);
      r(2) = r(2)*f(i);
   end
   i = i+1;
   j = j-1;
   k = k + 1;   
end
if i == j
   [m k] = min(r);
   r(k) = r(k) * f(i);
end
fo = [min(r) max(r)];   
m1 = fo(2) -fo(1);

% second way

for i = 1:numfact-1
    p(i,1) = prod(f(1:i));
    p(i,2) = prod(f(i+1:end));
end
[m2 k] = min(abs(p(:,1)-p(:,2)));




% pick the best one

if m1 < m2
    out = fo;
else
    out = [min(p(k,:)) max(p(k,:))];
end

