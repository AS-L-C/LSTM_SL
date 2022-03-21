function myAddText(ah,cx,cy,lloc,str,fs)
switch lloc(1)
    case 'b'
        ty = cy(1);
        val = 'bottom';
    case 't'
        ty = cy(2);
        val = 'top';
    case 'm'
        ty = mean(cy);
        val = 'middle';
end
switch lloc(2)
    case 'l'
        tx = cx(1);
        hal = 'left';
    case 'r'
        tx = cx(2);
        hal = 'right';
    case 'm'
        tx = mean(cx);
        hal = 'center';
end
text(ah,tx,ty,str,'FontSize',fs,'HorizontalAlignment',hal,'VerticalAlignment',val);
end