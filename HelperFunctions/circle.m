function h = circle(x,y,r,ax1,col)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(ax1,xunit, yunit,col);
end