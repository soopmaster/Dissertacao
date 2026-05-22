function draw_rectangle(Rec)
if size(Rec,1) == 2 && size(Rec,2) ==2
    rectangle('Position',[Rec(1,1),Rec(2,1),Rec(1,2)-Rec(1,1),Rec(2,2)-Rec(2,1)],'EdgeColor','k');
else
    disp('wrong dimention for rectangle');
end