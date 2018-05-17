function [htri] = plot_triangle (x, y, z, ttype, b)
% htri = triangle (x, y, ttype, b)
%   x, y are the coordinates of the node
%   ttpe = 1 is horizontal direction
%   ttpe = 2 is vertical direction
%   b is a blue tone color b=1,2,...

h = 1.0;
if b == 0, b = 1+rand; end
r = (1 - 1/b)/1.5;
g = (1 - 1/b)/1.0;
b = (1 + 1/b)/3.0;

htri_x = zeros(size(x));
htri_y = zeros(size(y));
htri_z = zeros(size(z));

if ttype == 1 % x
    x1 = x-(h*cos(pi/6));
    y1 = y+(h/2);
    y2 = y-(h/2);
    htri_x = arrayfun(@(x,x1,y,y1,z) ...
        line([x,x1],[y,y1],[z,z],'Color',[r g b],'LineWidth',2),x,x1,y,y1,z, 'UniformOutput', false);
    htri_y = arrayfun(@(x1,y1,y2,z) ...
        line([x1,x1],[y1,y2],[z,z],'Color',[r g b],'LineWidth',2),x1,y1,y2,z, 'UniformOutput', false);
    htri_z = arrayfun(@(x,x1,y,y2,z) ...
        line([x1,x],[y2,y],[z,z],'Color',[r g b],'LineWidth',2),x,x1,y,y2,z, 'UniformOutput', false);
elseif ttype == 2 % y
    z1 = z+(h/2);
    z2 = z-(h/2);
    y1 = y+(h*cos(pi/6));
    htri_x = arrayfun(@(x,y,y1,z,z1) ...
        line([x,x],[z,z1],[y,y1],'Color',[r g b],'LineWidth',2),x,y,y1,z,z1, 'UniformOutput', false);
    htri_y = arrayfun(@(x,y1,z1,z2) ...
        line([x,x],[z1,z2],[y1,y1],'Color',[r g b],'LineWidth',2),x,y1,z1,z2, 'UniformOutput', false);
    htri_z = arrayfun(@(x,y,y1,z,z2) ...
        line([x,x],[z2,z],[y1,y],'Color',[r g b],'LineWidth',2),x,y,y1,z,z2, 'UniformOutput', false);
elseif ttype == 3
    z1 = z-(h*cos(pi/6));
    y1 = y+(h/2);
    y2 = y-(h/2);
    htri_x = arrayfun(@(x,y,y1,z,z1) ...
        line([x,x],[z,z1],[y,y1],'Color',[r g b],'LineWidth',2),x,y,y1,z,z1, 'UniformOutput', false);
    htri_y = arrayfun(@(x,y1,y2,z1) ...
        line([x,x],[z1,z1],[y1,y2],'Color',[r g b],'LineWidth',2),x,y1,y2,z1, 'UniformOutput', false);
    htri_z = arrayfun(@(x,y,y2,z,z1) ...
        line([x,x],[z1,z],[y2,y],'Color',[r g b],'LineWidth',2),x,y,y2,z,z1, 'UniformOutput', false);
end

htri = [htri_x(:); htri_y(:); htri_z(:)];
