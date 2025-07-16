function width = fw50(x,y)
% only works for positive signals
N = length(y);
tot = sum(y);

[gar,centerindex]=max(y);

js = 1;

while ( (centerindex-js > 2) & (centerindex + js < N-2) & (sum(y(centerindex-js:centerindex+js))< tot/2) ) 
    js = js + 1;
end
if (centerindex+js<N-2) & (centerindex-js> 2)
width = x(centerindex+js)-x(centerindex-js);
else
    width = (x(end)-x(1))/2;

end

