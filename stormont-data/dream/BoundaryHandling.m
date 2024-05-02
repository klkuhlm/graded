function [x] = BoundaryHandling(x,ParRange,BoundHandling);
% Function to check whether parameter values remain within prior bounds

if strcmp(BoundHandling,'Reflect');
    [x] = ReflectBounds(x,ParRange);
end;
if strcmp(BoundHandling,'Bound');
    [x] = SetToBounds(x,ParRange);
end;
if strcmp(BoundHandling,'Fold');
    [x] = FoldBounds(x,ParRange);
end;