function [varargout]=plot_n_ellipsoid(center,P_inv,r,color,efficient)

% Edit of plot_ellipsoid function. Considers offset as well

% This function takes an ellipsoid P_inv matrix, a center, and a radius and plots the n-ellipsoid.  
% 
% If the 'efficient' flag is set, then what is passed to the function 
% via the second argument is the inv(sqrtm(P_inv)) instead of P_inv.    
% (JPK 12/2002)

% The shaded flag only works for 2 dimensional ellipsoids.
shaded = 0;

% Surface flag determines if the 3D plot is a line or a surface plot.
surface = 1;

if (nargin==5)&(efficient==1)
    A = P_inv;
elseif (nargin==4)|((nargin==5)&(efficient==0))
%     try
        A = inv(sqrtm(P_inv));
%     catch err
%         1
%     end
else
    error('Wrong number of input arguments or incorrect argument for efficiency flag.');
end

if size(P_inv,1)==2
   angle = linspace(0,2*pi);
   x_circ = r*cos(angle);
   y_circ = r*sin(angle);

   ellips = zeros(2,length(angle));
   for i=1:length(angle)
       ellips(:,i) = A*[ x_circ(i) ; y_circ(i) ] + center;
   end

   if shaded==1
      h = patch(ellips(1,:)',ellips(2,:)',color); 
      set(h,'EdgeColor',color)
   else
      h=plot(ellips(1,:),ellips(2,:),':','Color',color,'LineWidth',0.5);
       
   end
   
   if nargout>0
    varargout{1} = h;
   end
   
elseif size(P_inv,1)==3
   zangle = linspace(0,pi,50);
   angle = linspace(0,2*pi,50);
   x_sphere = [];
   y_sphere = [];
   z_sphere = [];

   if surface==1
      for i=1:length(zangle)
         x(i,:)=r*cos(angle)*sin(zangle(i));
         y(i,:)=r*sin(angle)*sin(zangle(i));
         z(i,1:length(angle))=r*cos(zangle(i));
      end

      for i=1:size(x,1)
         for j=1:size(x,2)
            transformed_point = A*[ x(i,j) ; y(i,j) ; z(i,j) ] + center;
            x_trans(i,j) = transformed_point(1);
            y_trans(i,j) = transformed_point(2);
            z_trans(i,j) = transformed_point(3);
         end
      end

      % surf(x_trans,y_trans,z_trans); 
      if nargout>0
          varargout{1} = mesh(x_trans,y_trans,z_trans);
      else
          mesh(x_trans,y_trans,z_trans,'FaceAlpha','0.5');
      end

      colormap([0 0 0]);
      % shading interp;
      % light;

   else
      for i=1:length(zangle)
         x_sphere(length(x_sphere)+1:length(x_sphere)+length(angle)) = r*cos(angle)*sin(zangle(i));
         y_sphere(length(y_sphere)+1:length(y_sphere)+length(angle)) = r*sin(angle)*sin(zangle(i));
         z_sphere(length(z_sphere)+1:length(z_sphere)+length(angle)) = r*cos(zangle(i));
      end

      ellips = zeros(3,length(x_sphere));
      for i=1:length(x_sphere)
         ellips(:,i) = A*[ x_sphere(i) ; y_sphere(i) ; z_sphere(i) ] + center;
      end

      if nargout>0
          varargout{1} = plot3(ellips(1,:),ellips(2,:),ellips(3,:));
      else
          plot3(ellips(1,:),ellips(2,:),ellips(3,:));
      end
       
   end

end

return
