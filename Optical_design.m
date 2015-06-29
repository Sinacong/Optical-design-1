function [image,focus,optics,probes,Mtot] = optical_design (beam,optics,probes,reflinv,plotflag)

%	function [image,focus,optics,probes,Mtot] = optical_design (beam,optics,probes,reflinv,plotflag)
%
%	general-purpose program to calculate the characteristics of an optical
%	  system in aid of optical design; calculations are done in the 
%	  first-order geometric optics limit but including third-order Seidel 
%	  aberrations; exact ray tracing is also performed for a number of
%	  standard rays
%	the system is assumed to be in vacuum/air (index of refraction = unity)
%	  and to possess planar symmetry (all optics are surfaces of rotation
%	  with axes in the x-y plane)
%	the optics implemented are lenses and reflectors: the thin-lens 
%	  approximation is used for lenses for the paraxial calculations, but
%	  the ray tracing is performed exactly
%	reflectors are defined by an index of refraction of -1 
%	off-axis optics are not treated separately in the paraxial calculations,
%	  i.e. they are considered to have the same properties as axial optics 
%	  once optical axis is re-oriented; a more precise treatment is not 
%	  available in general; the ray tracing on the other hand is exact
%	the user specifies the initial characteristics (at the "object") of a 
%	  light bundle and the characteristics and positions of an arbitrary 
%	  number of non-axial optics; all parameters are relative to the
%	  optical axis: when off-axis optics are present, the optical axis
%	  is propagated through them by ray tracing
%	four plots are generated:
%		1) axial and principal rays only, paraxial approximation
%		2) ray bundles from several points in the object, paraxial
%		   approximation
%		3) ray tracing for a set of 20 meridian and skew rays:
%			optical axis, axial marginal meridional, axial zonal
%			  meridional, principal meridional, edge marginal
%			  meridional (upper rim), edge marginal meridional 
%			  (lower rim), axial marginal sagittal, axial zonal
%			  sagittal, principal sagittal, edge marginal sagittal 
%			  (skew)
%			+ 10 similar rays with the object at infinity
%		4) ray tracing for a bundle of meridian and skew rays from 
%		   either:
%		      a single, user-selected point (not necessarily at the 
%			object); the bundle is made up of equally spaced 
%			tangential and sagittal rays
%		   or
%		      a multiple set of user-selected points; in this case the
%			bundle is made up of parallel (user-specified) rays
%	 	if the parameter reflinv is nonzero, an imaginary flat
%	  	  reflector is introduced after each real one to reverse the 
%		  first reflection, so that rays still propagate from left to 
%		  right (for ease of visualization)
%
%	input:
%	   beam = structure with fields:
%		w0 (m)		= half-width at object
%		theta0 (rad)	= half-divergence at object
%		stop (m)	= distance of aperture stop from object
%		aperture (m)	= half-width of aperture stop 
%		raysource_x,
%		raysource_y,
%		raysource_z	= point, specified with x (0 being the object),
%				  y and z coordinates as fractions of w0
%				  (0, 1/sqrt(2), 0 by default) from which
%				  ray bundle will be traced in plot 4
%				if these are vectors instead of scalars, the
%				  bundle is made up of parallel rays using the
%				  parameters raydir_y and raydir_z
%		raydir_y,
%		raydir_z	= y and z components (cosines) of unit vector
%				    defining ray propagation in the case of
%				    parallel rays (see above); 0 by default
%		   note: theta0 on one side and stop and aperture on the other
%			are mutually exclusive; if theta0 is specified, the 
%			aperture stop is the focal plane and its half-width is 
%			F*theta0; if theta0=[], stop and aperture are used
%	   optics (array) = structure with fields:
%		F (m)	= focal lengths of optical elements
%		d (m) 	= distance of each optical element from the previous one
%			  (for the first one, distance from object); distances
%			  are measured along the optical axis between second
%			  principal planes (defined by principal rays) (optic 
%			  surface in case of mirrors)
%		thickness (m) = thickness of optical element
%			for mirrors this is taken to be zero; for lenses, if 
%			  not specified it is defined to make the two surfaces
%			  meet at the edges of the ray bundle - this is also
%			  done if the thickness is specified but insufficient
%		C1 (m-1) = curvatures of entry surfaces (or sole surface in the
%			     case of a reflector)
%		C2 (m-1) = curvatures of exit surfaces
%		n 	= indexes of refraction
%		A2_1 (m-1) = for entry surfaces, coefficients in equation 
%		     x = C*R^2/(1+sqrt(1-C^2*R^2))+A2*R^2+A4*R^4 for aspherics
%		A2_2 (m-1) = A2 coefficients for exit surfaces
%		A4_1 (m-3) = for entry surfaces, coefficients in equation 
%		     x = C*R^2/(1+sqrt(1-C^2*R^2))+A2*R^2+A4*R^4 for aspherics
%		A4_2 (m-3) = A4 coefficients for exit surfaces
%		oad (m)	= off-axis distance: distance between entering optical 
%			    axis and axis of optic, perpendicularly to
%			    axis of optic, at intersection point of entering
%			    optical axis and first surface of optic (positive if
%			    optic is shifted upward)
%		oaa (rad) = off-axis angle: angle between entering optical axis
%			      and axis of optic (positive when axis of optic
%			      is shifted counterclockwise from initial axis)
%		shape	= type of optical element
%					0 = generic
%					1 = convex-plano lens
%					2 = spherical mirror
%					3 = paraboloidal mirror
%					4 = ellipsoidal mirror
%					5 = hyperboloidal mirror
%		   note: if shape=1 to 5 is specified, C2 is set to zero
%			if both F and n are given, C1 is ignored too
%		   note for lenses:
%			if F, C1 and n are provided, C2 is ignored (derived
%			  from F, C1 and thickness)
%			if F and n are provided and C1 and C2 are both empty,
%			  C2 is set to zero (convex-plano lens) and C1 derived 
%			  from F
%		   note for mirrors (n=-1):
%			C2 is set to zero
%			if F is provided, C1 is derived from it
%			note also that all the paraxial calculations are
%			  performed by using the approximate curvature and
%			  aspheric parameters; however, for the ray tracing
%			  the cases shape=3, 4, 5 are treated exactly
%				in the cases shape=4 and 5, A4_1 must be 
%				  provided; if it is not provided it is set to 
%				  zero, i.e. the spherical limit is assumed
%	   probes [array] = structure with fields:
%		x (m)	= vector of distances from object of additional
%			    locations at which some quantities are to be
%			    calculated
%	   reflinv [scalar] = if nonzero, an imaginary flat reflector is 
%				introduced after each real one to reverse the 
%		  		first reflection, so that rays still propagate 
%				from left to right (for ease of visualization)
%	   plotflag [scalar] = if zero, plots are not generated (default = 1)
%
%	output:
%	   image [array] = structure with fields
%		x (m)	= locations of images from each optical element,
%			    as distance from element (along the optical axis)
%		y (m)	= half-widths at images (located at image.x)
%		u (rad)	= angles of axial (marginal) rays 
%		    convention: positive when shifted counterclockwise from axis
%		up (rad)= angles of principal (chief) rays 
%		    convention: positive when shifted counterclockwise from axis
%		M	= magnification introduced by each optical element
%		aberrations = structure with fields; these are the total 
%		  aberrations at each image due to the combination of all the
%		  optics preceding it
%			TSC (m) = third-order transverse spherical aberration
%			SCC (m) = third-order sagittal coma
%			TCC (m) = third-order tangential coma
%			TPC (m) = third-order transverse Petzval curvature
%			TAC (m) = third-order transverse astigmatism
%			DC (m)  = third-order distortion
%			blur (m) = estimate of blur diameter, with spherical and
%			  coma aberrations minimized by image shift (note that
%			  since this includes the effects of field curvature and
%			  distortion, it will be larger than the blur directly
%			  visible by ray tracing)
%			OPDw (m) = optical path difference, worst-case scenario
%			  (no image shift)
%			OPDo (m) = optical path difference, with spherical and
%			  coma aberrations minimized by image shift
%	   focus [array] = structure with fields
%		x (m)	= locations of foci from each optical element,
%			    as distance from element (along the optical axis)
%		y (m)	= half-widths at foci
%		M	= magnification introduced by each optical element in
%			    propagating the focus
%		aberrations = structure with fields; these are the total 
%		  aberrations at each focus due to the combination of all the
%		  optics preceding it
%			TSC (m) = third-order transverse spherical aberration
%			SCC (m) = third-order sagittal coma
%			TCC (m) = third-order tangential coma
%			TAC (m) = third-order transverse astigmatism
%			TPC (m) = third-order transverse Petzval curvature
%			DC (m)  = third-order distortion
%			blur (m) = estimate of blur diameter, with spherical and
%			  coma aberrations minimized by image shift (note that
%			  since this includes the effects of field curvature and
%			  distortion, it will be larger than the blur directly
%			  visible by ray tracing)
%			OPDw (m) = optical path difference, worst-case scenario
%			  (no image shift)
%			OPDo (m) = optical path difference, with spherical and
%			  coma aberrations minimized by image shift
%	   optics [array] = structure with fields
%		h (m)	= half-widths at each optical element
%		u (rad)	= angles of axial (marginal) rays coming into optic
%		    convention: positive when shifted counterclockwise from axis
%		y (m) 	= heights of axial (marginal) rays
%		up (rad)= angles of principal (chief) rays coming into optic
%		    convention: positive when shifted counterclockwise from axis
%		yp (m) 	= heights of principal (chief) rays
%	   probes [array] = structure with fields
%		y (m) = half-widths at locations probes.x
%	   Mtot [scalar] = global magnification

%  Known problems:
%     infinite conjugates can have problems, particularly in the calculation of
%       aberrations
%     in practice it's best to detune the system imperceptibly to make the
%       infinites large but finite
%

% input checks and variable packing
input_variables={'beam','optics','probes'};
l_in=length(input_variables);

% fields in alphabetical order (essential here)
% note that structure of the program requires that field names be all distinct
%   (even across variables)
input_fields{1}={'aperture','raydir_y','raydir_z','raysource_x', ...
		 'raysource_y','raysource_z','stop','theta0','w0'};
input_fields{2}={'A2_1','A2_2','A4_1','A4_2','C1', ...
		 'C2','F','d','n','oaa','oad','shape','thickness'};
input_fields{3}={'x'};

for i=l_in:-1:1,
% check for empty or non-structure input variables and correct
 eval(['if nargin<',int2str(i),'|isempty(',input_variables{i}, ...
          ')|~isstruct(',input_variables{i},'),', ...
  		input_variables{i},'=struct(input_fields{', ...
			int2str(i),'}{1},NaN);', ...
       'end']); 
% transfer input variables to cell column array
 eval(['input_data{',int2str(i),'}=',input_variables{i},';']);
 input_data{i}=input_data{i}(:);
% check for necessary fields that have not been defined
 input_filled{i}=isfield(input_data{i},input_fields{i});
% remove extraneous fields
 for fi=fieldnames(input_data{i})'
  if ~strcmp(fi,input_fields{i}),
   input_data{i}=rmfield(input_data{i},fi);
  end
 end
end

for i=1:l_in,
 for j=1:length(input_filled{i}),
  if ~input_filled{i}(j),
% create fields that have not been defined yet
   for k=1:length(input_data{i}),
    input_data{i}=setfield(input_data{i},{k},input_fields{i}{j},NaN);
   end
  end
% empty fields cause problems, so we turn them to NaNs
  for k=1:length(input_data{i}),
   if isempty(getfield(input_data{i},{k},input_fields{i}{j})),
    input_data{i}=setfield(input_data{i},{k},input_fields{i}{j},NaN);
   end
  end
 end
end

% unfold parameters
for i=1:l_in,
% transfer inputs to cell array of cell matrices
  cell_data{i}=struct2cell(orderfields(input_data{i}));
% transfer to individual variables named after the fields
  for j=1:length(input_fields{i}),
    eval([input_fields{i}{j},'=[cell_data{i}{j,:}];']);
  end
end
% defaults and incompatible input parameters
if ~exist('reflinv')|isempty(reflinv)|isnan(reflinv),reflinv=0;end
if ~exist('plotflag')|isempty(plotflag)|isnan(plotflag),plotflag=1;end
if isnan(w0),
  w0=0.035;
  disp(['default half-width ',num2str(w0),' m']);
end
if isnan(theta0)&(isnan(stop)|isnan(aperture)),
  theta0=0.01;
  disp(['default half-divergence ',num2str(theta0),' rad']);
end
[x,probesorder]=sort(x);
[dummy,probesorder]=sort(probesorder);
if isnan(raysource_x),
  raysource_x=0;
end
if isnan(raysource_y),
  raysource_y=1/sqrt(2);
end
if isnan(raysource_z),
  raysource_z=0;
end
if isnan(raydir_y),
  raydir_y=0;
end
if isnan(raydir_z),
  raydir_z=0;
end
raydir_y=raydir_y(1);
raydir_z=raydir_z(1);
if raydir_y^2+raydir_z^2>1-5*eps,
  raydir_z=0;
end
raydir_x=sqrt(1-raydir_y^2-raydir_z^2);
mray4=min([length(raysource_x),length(raysource_y),length(raysource_z)]);
raysource_x=raysource_x(1:mray4);
raysource_y=raysource_y(1:mray4);
raysource_z=raysource_z(1:mray4);
raysource_x=raysource_x(:)';
raysource_y=raysource_y(:)';
raysource_z=raysource_z(:)';
shape(isnan(shape))=0;
thickness(isnan(thickness))=0;
% default: on-axis optics
oaa(isnan(oaa))=0;
oad(isnan(oad))=0;
% default index of refraction = 1 (air lenses)
n(isnan(n))=1+eps;
A2_1(isnan(A2_1))=0;
A2_2(isnan(A2_2))=0;
A4_1(isnan(A4_1))=0;
A4_2(isnan(A4_2))=0;
t1=find(shape>1&shape<6);
n(t1)=-1;
thickness(t1)=0;
t1=find(shape>0&shape<6&~isnan(F));
C1(t1)=1./F(t1)./(n(t1)-1);
% convex-plano lenses
t1=find(shape==1);
C2(t1)=0;A2_1(t1)=0;A2_2(t1)=0;A4_1(t1)=0;A4_2(t1)=0;
% spherical mirrors
t1=find(shape==2);
n(t1)=-1;A2_1(t1)=0;A2_2(t1)=0;A4_1(t1)=0;A4_2(t1)=0;
% paraboloidal mirrors
t1=find(shape==3);
n(t1)=-1;A2_1(t1)=0;A2_2(t1)=0;A4_1(t1)=-C1(t1).^3/8;A4_2(t1)=0;
% ellipsoidal or hyperboloidal mirrors
t1=find(shape==4|shape==5);
n(t1)=-1;A2_1(t1)=0;A2_2(t1)=0;A4_2(t1)=0;
% lenses
t1=find(n~=-1&~isnan(F)&~isnan(C1));
C2(t1)=(-1./F(t1)./(n(t1)-1)+C1(t1)+2*A2_1(t1))./ ...
       (1-thickness(t1).*(C1(t1)+2*A2_1(t1)).*(1-1./n(t1)))-2*A2_2(t1);
t1=find(n~=-1&~isnan(F)&isnan(C1)&isnan(C2));
C2(t1)=0;
C1(t1)=(1./F(t1)./(n(t1)-1)+2*A2_2(t1))./ ...
       (1+thickness(t1)*2.*A2_2(t1).*(1-1./n(t1)))-2*A2_1(t1);
% mirrors
t1=find(n==-1);
C2(t1)=0;
t1=find(n==-1&~isnan(F));
C1(t1)=1./F(t1)./(n(t1)-1)-2*(A2_1(t1)-A2_2(t1));

% calculate focal lengths
F=1./(n-1)./(C1-C2+2*(A2_1-A2_2)+thickness.*(C1+2*A2_1).*(C2+2*A2_2).*(1-1./n));

dist=cumsum(d);

% First pass at paraxial calculations just to determine approximate y extent of
%   optics

% Calculate focal lengths projected on new optical axes (for off-axis optics)
Fr=sqrt(F.^2.+(F.*tan(oaa)-oad).^2);

x_image(1)=1/(1/Fr(1)-1/d(1));
% for theta and y_image we do not take absolute values here because they
%   have to remain compatible for later sums and differences (upon imaging,
%   up-down symmetric bundle becomes asymmetric, so an initial +-theta0 becomes
%   theta1 and theta2; note however that when inverting both the initial w0 and
%   theta, the resulting theta1, theta2, and y_image are all inverted, at all
%   successive stages); we'll take the absolute value at the end
% in the case of stops and apertures, everything is symmetric and so we take
%   absolute values as we go along
if ~isnan(theta0),
  theta(:,1)=[1;-1]*theta0(1)*(1-d(1)/Fr(1))-w0(1)/Fr(1);
  h_optics(1)=w0(1)+theta0(1)*d(1);
  if n(1)==-1&~reflinv,
    theta(:,1)=-theta(:,1);
  end
else
% stops and apertures are the successive images of the initial stop and aperture
  stops(1)=1/(1/Fr(1)-1/(d(1)-stop(1)));
  apertures(1)=abs(aperture(1)*Fr(1)/(d(1)-stop(1)-Fr(1)));
  h_optics(1)=max(abs(aperture(1)*d(1)/stop(1)+[1;-1]*w0(1)*(d(1)/stop(1)-1)));
end
y_image(1)=-w0(1)*Fr(1)/(d(1)-Fr(1));
if n(1)==-1&~reflinv,
  y_image(1)=-y_image(1);
end

if length(d)~=1,
  if ~isnan(theta0),
    h_optics(2)=w0(1)*abs(1-d(2)/Fr(1))+theta0(1)*abs(dist(2)-d(1)*d(2)/Fr(1));
  else
    h_optics(2)=max(abs(repmat(apertures(1)*(d(2)-x_image(1))/ ...
				(stops(1)-x_image(1)),2,1)+ ...
	     [1;-1]*y_image(1)*((d(2)-x_image(1))/(stops(1)-x_image(1))-1)));
  end
end
for i=2:length(d),
  x_image(i)=1/(1/Fr(i)-1/(d(i)-x_image(i-1)));
  if ~isnan(theta0),
    theta(:,i)=theta(:,i-1)*(1-(d(i)-x_image(i-1))/Fr(i))-y_image(i-1)/Fr(i);
    if n(i)==-1&~reflinv,
      theta(:,i)=-theta(:,i);
    end
  else
    stops(i)=1/(1/Fr(i)-1/(d(i)-stops(i-1)));
    apertures(i)=abs(apertures(i-1)*Fr(i)/(d(i)-stops(i-1)-Fr(i)));
  end
  y_image(i)=-y_image(i-1)*Fr(i)/(d(i)-x_image(i-1)-Fr(i));
  if n(i)==-1&~reflinv,
    y_image(i)=-y_image(i);
  end
  if i~=length(d),
    if ~isnan(theta0),
      h_optics(i+1)=max(abs(y_image(i-1)*(1-d(i+1)/Fr(i))+ ...
     		theta(:,i-1)*(d(i+1)+d(i)-x_image(i-1)-d(i+1)* ...
			(d(i)-x_image(i-1))/Fr(i))));
    else
      h_optics(i+1)=max(abs(repmat(apertures(i)*(d(i+1)-x_image(i))/ ...
				      (stops(i)-x_image(i)),2,1)+ ...
		[1;-1]*y_image(i)*((d(i+1)-x_image(i))/ ...
				   (stops(i)-x_image(i))-1)));
    end
  end
end


% initialize rays for ray tracing
% rays = 3 x rays x surfaces (x,y,z coordinates)
% raydir = 3 x rays x surfaces (x,y,z direction cosines)
% the first ray must always be the optical axis
if mray4==1,
  nray4=7;
  if ~isnan(theta0),
    arr4=linspace(-1,1,nray4)*(nray4-1)/2;
    raydir_x=[cos(theta0(1)*arr4/nray4),cos(theta0(1)*arr4/nray4)];
    raydir_y=[sin(theta0(1)*arr4/nray4),zeros(1,nray4)];
    raydir_z=[zeros(1,nray4),sin(theta0(1)*arr4/nray4)];
  else
    thetasource_y=linspace(atan((-aperture(1)-w0(1)*raysource_y)/ ...
			      (stop(1)-raysource_x)), ...
			 atan((aperture(1)-w0(1)*raysource_y)/ ...
			      (stop(1)-raysource_x)),nray4);
    thetasource_z=linspace(atan((-aperture(1)-w0(1)*raysource_z)/ ...
			      (stop(1)-raysource_x)), ...
			 atan((aperture(1)-w0(1)*raysource_z)/ ...
			      (stop(1)-raysource_x)),nray4);
    raydir_x=cos([thetasource_y,thetasource_z]);
    raydir_y=[sin(thetasource_y),zeros(1,nray4)];
    raydir_z=[zeros(1,nray4),sin(thetasource_y)];
  end
end
if ~isnan(theta0),
% with the object at infinity, the stop cannot be at infinity as well, so for 
%   this application we place it at the real object
  rays=[[0;0;0], ...		%optical axis
        [0;0;0], ...		%axial marginal meridional
        [0;0;0], ...		%axial zonal meridional
        [0;w0(1);0], ...	%principal meridional
        [0;w0(1);0], ...	%edge marginal meridional (upper rim)
        [0;w0(1);0], ...	%edge marginal meridional (lower rim)
        [0;0;0], ...		%axial marginal sagittal
        [0;0;0], ...		%axial zonal sagittal
        [0;0;w0(1)], ...	%principal sagittal
        [0;w0(1);0], ...	%edge marginal sagittal (skew)
	[0;0;0], ...	%same for object at infinity
        [0;w0(1);0], ...
        [0;w0(1)/sqrt(2);0], ...
        [0;0;0], ...
        [0;w0(1);0], ...
        [0;-w0(1);0], ...
        [0;0;w0(1)], ...
        [0;0;w0(1)/sqrt(2)], ...
        [0;0;0], ...
        [0;w0(1);0]];
  nray3=size(rays,2);
  raydir=[[1;0;0], ...
	    [cos(theta0(1));sin(theta0(1));0], ...
	    [1/sqrt(1+(tan(theta0(1)))^2/2); ...
	     1/sqrt(1+2/(tan(theta0(1)))^2);0], ...
	    [1;0;0], ...
	    [cos(theta0(1));sin(theta0(1));0], ...
	    [cos(theta0(1));-sin(theta0(1));0], ...
	    [cos(theta0(1));0;sin(theta0(1))], ...
	    [1/sqrt(1+(tan(theta0(1)))^2/2);0; ...
	     1/sqrt(1+2/(tan(theta0(1)))^2)], ...
	    [1;0;0], ...
	    [cos(theta0(1));0;sin(theta0(1))], ...
	    [1;0;0], ...
	    [1;0;0], ...
	    [1;0;0], ...
	    [cos(theta0(1));sin(theta0(1));0], ...
	    [cos(theta0(1));sin(theta0(1));0], ...
	    [cos(theta0(1));sin(theta0(1));0], ...
	    [1;0;0], ...
	    [1;0;0], ...
	    [cos(theta0(1));0;sin(theta0(1))], ...
	    [cos(theta0(1));sin(theta0(1));0]];
else
% with the object at infinity, the principal rays are ill-defined, so we trace 
%   them from the (real) object edge
  rays=[[0;0;0], ...		%optical axis
        [0;0;0], ...		%axial marginal meridional
        [0;0;0], ...		%axial zonal meridional
        [0;w0(1);0], ...	%principal meridional
        [0;w0(1);0], ...	%edge marginal meridional (upper rim)
        [0;w0(1);0], ...	%edge marginal meridional (lower rim)
        [0;0;0], ...		%axial marginal sagittal
        [0;0;0], ...		%axial zonal sagittal
        [0;0;w0(1)], ...	%principal sagittal
        [0;w0(1);0], ...	%edge marginal sagittal (skew)
	[0;0;0], ...	%same for object at infinity
        [0;aperture(1);0], ...
        [0;aperture(1)/sqrt(2);0], ...
        [0;w0(1);0], ...
        [0;w0(1);0], ...
        [0;w0(1);0], ...
        [0;0;aperture(1)], ...
        [0;0;aperture(1)/sqrt(2)], ...
        [0;0;w0(1)], ...
        [0;w0(1);0]];
  nray3=size(rays,2);
  raydir=[[1;0;0], ...
	    [stop(1);aperture(1);0], ...
	    [stop(1);aperture(1)/sqrt(2);0], ...
	    [stop(1);-w0(1);0], ...
	    [stop(1);aperture(1)-w0(1);0], ...
	    [stop(1);-aperture(1)-w0(1);0], ...
	    [stop(1);0;aperture(1)], ...
	    [stop(1);0;aperture(1)/sqrt(2)], ...
	    [stop(1);0;-w0(1)], ...
	    [stop(1);-w0(1);aperture(1)], ...
	    [1;0;0], ...
	    [1;0;0], ...
	    [1;0;0], ...
	    [stop(1);-w0(1);0], ...
	    [stop(1);aperture(1)-w0(1);0], ...
	    [stop(1);-aperture(1)-w0(1);0], ...
	    [1;0;0], ...
	    [1;0;0], ...
	    [stop(1);0;-w0(1)], ...
	    [stop(1);-w0(1);aperture(1)]];
  raydir=raydir./repmat(sqrt(sum(raydir.^2)),3,1);
end
% rays for plot 4
if mray4>1,
  if raydir_x>5*eps,
    l=-w0(1)*raysource_x/raydir_x;
  else
    l=0;
  end
  rays=[rays,w0(1)*[raysource_x;raysource_y;raysource_z]+ ...
	       [raydir_x;raydir_y;raydir_z]*l];
  raydir=[raydir,repmat([raydir_x;raydir_y;raydir_z],1,mray4)];
else
  l=-w0(1)*raysource_x./raydir_x;
  rays=[rays,w0(1)*repmat([raysource_x;raysource_y;raysource_z], ...
		1,2*nray4)+[raydir_x;raydir_y;raydir_z].*repmat(l,3,1)];
  raydir=[raydir,[raydir_x;raydir_y;raydir_z]];
end
if isnan(theta0),
  y_ray_stop=rays(2,:)+(stop(1)-rays(1,:))*raydir(2,:)./raydir(1,:);
end
% make rays travel to the right
for i=1:size(raydir,2),
  if raydir(1,i)<0,raydir(:,i)=-raydir(:,i);end
end
raytracing(:,:,1)=rays;

% ray tracing
opt_axis(1)=0;
for i=1:length(d),
  pri1(i)=0;
% first surface

% move to coordinate system with origin in surface vertex and x axis along
%   optic's symmetry axis
  if shape(i)>0&shape(i)<6,
    alpha_opt(1,i)=C1(i)/2;
    alpha_opt(2,i)=C2(i)/2;
    if abs(alpha_opt(1,i))<5*eps,
      beta_opt(1,i)=0;
    else
      beta_opt(1,i)=A4_1(i)/alpha_opt(1,i)^2+alpha_opt(1,i);
    end
    if abs(alpha_opt(2,i))<5*eps,
      beta_opt(2,i)=0;
    else
      beta_opt(2,i)=A4_2(i)/alpha_opt(2,i)^2+alpha_opt(2,i);
    end
% xshift is the distance along the surface axis between the vertex and the
%   projection of the intersection of the original optical axis with the surface
    if abs(beta_opt(1,i))<5*eps,
      xshift(i)=alpha_opt(1,i)*oad(i)^2;
    else
      discrim=1-4*alpha_opt(1,i)*beta_opt(1,i)*oad(i)^2;
      if discrim<-5*eps,
	disp(['No intersection with first surface of optic ',int2str(i)]);
	return
      elseif discrim<0,
	discrim=0;
      end
      xshift(i)=0.5/beta_opt(1,i)*(1-sqrt(discrim));
    end
% second-degree equation to find upper end of optic
    acoe=alpha_opt(1,i)*(tan(oaa(i)))^2+beta_opt(1,i);
    bcoe=2*alpha_opt(1,i)*tan(oaa(i))* ...
	 (oad(i)-xshift(i)*tan(oaa(i))-h_optics(i)/cos(oaa(i)))-1;
    ccoe=alpha_opt(1,i)* ...
	 (oad(i)-xshift(i)*tan(oaa(i))-h_optics(i)/cos(oaa(i)))^2;
    if acoe,
      discrim=bcoe^2-4*acoe*ccoe;
      if discrim<-5*eps,
	disp(['No intersection with first surface of optic ',int2str(i)]);
	return
      elseif discrim<0,
	discrim=0;
      end
      xup=0.5*(-bcoe+sqrt(discrim)*[1;-1])/acoe;
      [dummy,jup]=min(abs(xup-xshift(i)));
      xup=xup(jup);
    else
      xup=-ccoe/bcoe;
    end
% second-degree equation to find lower end of optic
    acoe=alpha_opt(1,i)*(tan(oaa(i)))^2+beta_opt(1,i);
    bcoe=2*alpha_opt(1,i)*tan(oaa(i))* ...
	 (oad(i)-xshift(i)*tan(oaa(i))+h_optics(i)/cos(oaa(i)))-1;
    ccoe=alpha_opt(1,i)* ...
	 (oad(i)-xshift(i)*tan(oaa(i))+h_optics(i)/cos(oaa(i)))^2;
    if acoe,
      discrim=bcoe^2-4*acoe*ccoe;
      if discrim<-5*eps,
	disp(['No intersection with first surface of optic ',int2str(i)]);
	return
      elseif discrim<0,
	discrim=0;
      end
      xdo=0.5*(-bcoe+sqrt(discrim)*[1;-1])/acoe;
      [dummy,jdo]=min(abs(xdo-xshift(i)));
      xdo=xdo(jdo);
    else
      xdo=-ccoe/bcoe;
    end
% left and right points
    xleft=min(xup,xdo);
    xright=max(xup,xdo);
  else
    discrim=1-(C1(i)*oad(i))^2;
    if discrim<-5*eps,
      disp(['No intersection with first surface of optic ',int2str(i)]);
      return
    elseif discrim<0,
      discrim=0;
    end
    xshift(i)=C1(i)*oad(i)^2/(1+sqrt(discrim))+ ...
	      A2_1(i)*oad(i)^2+A4_1(i)*oad(i)^4;
    xup=fzero(@(x) ...
	aspheric_skew_intersects(x,h_optics(i),xshift(i),oad(i),oaa(i), ...
				 C1(i),A2_1(i),A4_1(i)),xshift(i));
    xdo=fzero(@(x) ...
	aspheric_skew_intersects(x,-h_optics(i),xshift(i),oad(i),oaa(i), ...
				 C1(i),A2_1(i),A4_1(i)),xshift(i));
    xleft=min(xup,xdo);
    xright=max(xup,xdo);
  end
  ytemp=(-[xup;xdo]+xshift(i))*tan(oaa(i))+[1;-1]*h_optics(i)/cos(oaa(i)) ...
	-oad(i);
  y_adjust(i)=mean(ytemp)+oad(i);
  h_optics(i)=abs(diff(ytemp))/2;
  if n(i)~=-1,
% second surface
    ytemp1=y_adjust(i)+h_optics(i);
    ytemp2=y_adjust(i)-h_optics(i);
    if shape(i)>0&shape(i)<6,
      if beta_opt(2,i),
        discrim1=1-4*alpha_opt(2,i)*beta_opt(2,i)*ytemp1^2;
        discrim2=1-4*alpha_opt(2,i)*beta_opt(2,i)*ytemp2^2;
        if discrim1<-5*eps,
	  disp(['No intersection with second surface of optic ',int2str(i)]);
	  return
        elseif discrim1<0,
	  discrim1=0;
        end
        if discrim2<-5*eps,
	  disp(['No intersection with second surface of optic ',int2str(i)]);
	  return
        elseif discrim2<0,
	  discrim2=0;
        end
        xtemp1=0.5/beta_opt(2,i)*(1+[1;-1]*sqrt(discrim1));
        [dummy,jx]=min(abs(xtemp1));
        xtemp1=xtemp1(jx);
	xtemp2=0.5/beta_opt(2,i)*(1+[1;-1]*sqrt(discrim2));
        [dummy,jx]=min(abs(xtemp2));
        xtemp2=xtemp2(jx);
      else
	xtemp1=alpha_opt(2,i)*ytemp1^2;
	xtemp2=alpha_opt(2,i)*ytemp2^2;
      end
    else
      discrim1=1-(C2(i)*ytemp1)^2;
      discrim2=1-(C2(i)*ytemp2)^2;
      if discrim1<-5*eps,
	disp(['No intersection with second surface of optic ',int2str(i)]);
	return
      elseif discrim1<0,
        discrim1=0;
      end
      if discrim2<-5*eps,
	disp(['No intersection with second surface of optic ',int2str(i)]);
	return
      elseif discrim2<0,
	discrim2=0;
      end
      xtemp1=C2(i)*ytemp1^2/(1+sqrt(discrim1))+ ...
  	     A2_2(i)*ytemp1^2+A4_2(i)*ytemp1^4;
      xtemp2=C2(i)*ytemp2^2/(1+sqrt(discrim2))+ ...
	     A2_2(i)*ytemp2^2+A4_2(i)*ytemp2^4;
    end
    if C1(i)*C2(i)>=0,
      if abs(C1(i))>abs(C2(i)),
        minthick=xright-max(xtemp1,xtemp2);
      else
        minthick=xleft-min(xtemp1,xtemp2);
      end
    else
      minthick=xright-min(xtemp1,xtemp2);
    end
    thickness(i)=max(thickness(i),minthick);
  end
% distance from first surface vertex to second principal plane (calculated using
%   aspheric expansion)
  pri1(i)=thickness(i)*(1-(C1(i)+2*A2_1(i))/n(i)/ ...
       (C1(i)+2*A2_1(i)-C2(i)-2*A2_2(i)+thickness(i)* ...
		(C1(i)+2*A2_1(i))*(C2(i)+2*A2_2(i))*(1-1/n(i))));
% in the case of a mirror the reference point is the impact point of the 
%   optical axis
  if n(i)==-1,
    pri1(i)=xshift(i);
  end
  rays_con=[cos(oaa(i)),sin(oaa(i)),0; ...
	    -sin(oaa(i)),cos(oaa(i)),0; ...
	    0,0,1] * (rays(:,:,i)-repmat( ...
	[dist(i)-(pri1(i)-xshift(i))/cos(oaa(i));0;0],1,size(rays,2))) ...
		     + repmat([xshift(i);-oad(i);0],1,size(rays,2));
  raydir_con=[cos(oaa(i)),sin(oaa(i)),0; ...
	      -sin(oaa(i)),cos(oaa(i)),0; ...
	      0,0,1] * raydir(:,:,i);
% store coordinates of initial system origin
  origin(:,i)=[cos(oaa(i)),sin(oaa(i)),0; ...
	       -sin(oaa(i)),cos(oaa(i)),0; ...
	       0,0,1] * [-dist(i)+(pri1(i)-xshift(i))/cos(oaa(i));0;0] ...
	 + [xshift(i);-oad(i);0];
% opt_axis(2*i) is the cumulated distance along the optical axis up to the
%   second principal point of this optic
  dummy=(pri1(i)-rays_con(1,1))/raydir_con(1,1);
  opt_axis(2*i)=norm(origin(:,i)-rays_con(:,1)-raydir_con(:,1)*dummy);
% store direction vector of initial optical axis
  opt_axis_ray=raydir_con(:,1);
% get intersection points and surface normals
  [l,s_nor] = ray_aspheric_intersect( ...
		rays_con,raydir_con,shape(i),C1(i),A2_1(i),A4_1(i),d(i));
  raytracing(:,:,2*i)=[cos(oaa(i)),-sin(oaa(i)),0; ...
		       sin(oaa(i)),cos(oaa(i)),0;
		       0,0,1]*(rays_con+raydir_con.*repmat(l,3,1)+ ...
     repmat([-xshift(i);oad(i);0],1,size(rays_con,2)))+ ...
     repmat([dist(i)-(pri1(i)-xshift(i))/cos(oaa(i));0;0],1,size(rays_con,2));
% angles of incidence
  incc=dot(raydir_con,s_nor);
  rays_con=rays_con+raydir_con.*repmat(l,3,1);
  if incc,
    raydir_con=s_nor.*repmat(incc.*(sqrt((n(i)^2-1)./incc.^2+1)/n(i) ...
			-1/abs(n(i))),3,1)+raydir_con/abs(n(i));
  else
    raydir_con=s_nor/n(i).*repmat(sqrt(n(i)^2-1),3,1)+raydir_con/abs(n(i));
  end
% mirror: special case
  if n(i)==-1,
    if reflinv,
% mirror: flip everything from left to right
      rays_con(1,:)=2*xshift(i)-rays_con(1,:);
      raydir_con(1,:)=-raydir_con(1,:);
   end
   raytracing(:,:,2*i+1)=rays_con;
  else
% move to coordinate system of second surface
    rays_con=rays_con-repmat([thickness(i);0;0],1,size(rays_con,2));
    origin(1,i)=origin(1,i)-thickness(i);
% get intersection points and surface normals
    [l,s_nor] = ray_aspheric_intersect( ...
         rays_con,raydir_con,shape(i),C2(i),A2_2(i),A4_2(i), ...
	 -xshift(i)+thickness(i)+(xtemp2+xtemp1)/2);
% ray points in the frame of reference of the second surface
    raytracing(:,:,2*i+1)=rays_con+raydir_con.*repmat(l,3,1);
% angles of incidence
    incc=dot(raydir_con,s_nor);
    rays_con=rays_con+raydir_con.*repmat(l,3,1);
    if incc,
      raydir_con=s_nor.*repmat(incc.*(sqrt((1/n(i)^2-1)./ ...
			incc.^2+1)*n(i)-abs(n(i))),3,1)+raydir_con*abs(n(i));
    else
      raydir_con=s_nor*n(i).*repmat(sqrt(1/n(i)^2-1),3,1)+raydir_con*abs(n(i));
    end
  end
% distance from old origin to second principal plane along entering optical 
%  axis, plus distance from second principal plane to intersection with second 
%   surface along outgoing optical axis
  if abs(raydir_con(1,1))>5*eps,
    opt_axis_length=opt_axis(2*i)+ ...
		  (thickness(i)-pri1(i)+rays_con(1,1))/raydir_con(1,1);
  else
    opt_axis_length=opt_axis(2*i);
  end
% rotation angle from current system to new optical axis
  rotan=atan2(raydir_con(2,1),raydir_con(1,1));
% move to coordinate system whose x-axis is the new optical axis; the origin
%   is the initial origin as seen in this new system, with the following
%   convention: the distance from the origin to a point along the optical
%   axis is equal to the sum of the optical path segment lengths from second
%   principal point to second principal point (note that the segments may be 
%   discontinuous)
% we put the second surface ray tracing points in the new coordinate system, as
%   well as the outcoming rays
  raytracing(:,:,2*i+1)=[cos(rotan),sin(rotan),0; ...
			 -sin(rotan),cos(rotan),0; ...
			 0,0,1]* ...
	(raytracing(:,:,2*i+1)-repmat(rays_con(:,1),1,size(raytracing,2)))+ ...
	repmat([opt_axis_length;0;0],1,size(raytracing,2));
  rays(:,:,i+1)=[cos(rotan),sin(rotan),0; ...
		 -sin(rotan),cos(rotan),0; ...
		 0,0,1]* ...
	(rays_con-repmat(rays_con(:,1),1,size(rays_con,2)))+ ...
	repmat([opt_axis_length;0;0],1,size(raytracing,2));
  raydir(:,:,i+1)=[cos(rotan),sin(rotan),0; ...
		   -sin(rotan),cos(rotan),0;
		   0,0,1]*raydir_con;
% new origin as seen from old origin
  origin(:,i)=[cos(rotan),sin(rotan),0; ...
	       -sin(rotan),cos(rotan),0; ...
	       0,0,1]*(rays_con(:,1)-origin(:,i))-[opt_axis_length;0;0];
  pivot_an(i)=rotan-atan2(opt_axis_ray(2),opt_axis_ray(1));
  opt_axis(2*i+1)=opt_axis(2*i);
end


% calculate focal lengths projected on new optical axes (for off-axis optics);
%   to obtain as good as possible a match with the ray tracing, we use a
%   higher-order formula here and not just the paraxial one
im=find(n==-1);
it=find(n~=-1);
Fr(it)=sqrt(F(it).^2.+(F(it).*tan(oaa(it)) ...
	    -oad(it)-(pri1(it)-xshift(it)).*tan(oaa(it))).^2);
Fr(im)=sqrt((F(im)+xshift(im)).^2.+ ...
	    ((F(im)+xshift(im)).*tan(oaa(im))-oad(im)).^2);

x_image(1)=1/(1/Fr(1)-1/d(1));
x_focus(1)=Fr(1);
imax=max(find(x<=dist(1)));
if isempty(imax),imax=0;end
% for theta and y_image we do not take absolute values here because they
%   have to remain compatible for later sums and differences (upon imaging,
%   up-down symmetric bundle becomes asymmetric, so an initial +-theta0 becomes
%   theta1 and theta2; note however that when inverting both the initial w0 and
%   theta, the resulting theta1, theta2, and y_image are all inverted, at all
%   successive stages); we'll take the absolute value at the end
% in the case of stops and apertures, everything is symmetric and so we take
%   absolute values as we go along
if ~isnan(theta0),
  theta(:,1)=[1;-1]*theta0(1)*(1-d(1)/Fr(1))-w0(1)/Fr(1);
  if n(1)==-1&~reflinv,
    theta(:,1)=-theta(:,1);
  end
  y_focus(1)=abs(theta0(1)*abs(Fr(1)));
  h_optics(1)=w0(1)+theta0(1)*d(1);
  y_probes(1:imax)=w0(1)+theta0(1)*x(1:imax);
% axial and principal ray parameters for object, coming into optic
  ray_opt(:,1,1)=[theta0(1);theta0(1)*d(1)];
  ray_opt(:,2,1)=[0;w0(1)];
% axial and principal ray parameters for an imaginary object at infinity; the
%   stop then cannot be at infinity as well, so for this application we place it
%   at the real object
  ray_opt(:,3,1)=[0;w0(1)];
  ray_opt(:,4,1)=[theta0(1);0];
else
% stops and apertures are the successive images of the initial stop and aperture
  stops(1)=1/(1/Fr(1)-1/(d(1)-stop(1)));
  apertures(1)=abs(aperture(1)*Fr(1)/(d(1)-stop(1)-Fr(1)));
  y_focus(1)=(abs(aperture(1))+abs(w0(1)))*abs(Fr(1)/stop(1));
  h_optics(1)=max(abs(aperture(1)*d(1)/stop(1)+[1;-1]*w0(1)*(d(1)/stop(1)-1)));
  y_probes(1:imax)=max(abs(aperture(1)*repmat(x(1:imax),2,1)/stop(1)+ ...
		    [1;-1]*w0(1)*(x(1:imax)/stop(1)-1)));
% axial and principal ray parameters for object
  ray_opt(:,1,1)=[aperture(1)/stop(1);aperture(1)/stop(1)*d(1)];
  ray_opt(:,2,1)=[-w0(1)/stop(1);w0(1)*(1-d(1)/stop(1))];
% axial and principal ray parameters for an imaginary object at infinity; the
%   principal ray is ill-defined in this case, so we trace it from the object
%   edge
  ray_opt(:,3,1)=[0;aperture(1)];
  ray_opt(:,4,1)=[-w0(1)/stop(1);w0(1)*(1-d(1)/stop(1))];
end
y_image(1)=-w0(1)*Fr(1)/(d(1)-Fr(1));
if n(1)==-1&~reflinv,
  y_image(1)=-y_image(1);
end

% axial and principal ray parameters for image
ray_im(:,1:2,1)=[1,0;x_image(1),1]*[1,-1/Fr(1);0,1]*ray_opt(:,1:2,1);
ray_foc(:,1:2,1)=[1,0;x_focus(1),1]*[1,-1/Fr(1);0,1]*ray_opt(:,3:4,1);
if n(1)==-1&~reflinv,
  ray_im(:,1:2,1)=-ray_im(:,1:2,1);
  ray_foc(:,1:2,1)=-ray_foc(:,1:2,1);
end
imin=imax+1;
if length(d)==1,
  imax=length(x);
else
  imax=max(find(x<=dist(2)));
  if isempty(imax),imax=0;end
  if ~isnan(theta0),
    h_optics(2)=w0(1)*abs(1-d(2)/Fr(1))+theta0(1)*abs(dist(2)-d(1)*d(2)/Fr(1));
  else
    h_optics(2)=max(abs(repmat(apertures(1)*(d(2)-x_image(1))/ ...
				(stops(1)-x_image(1)),2,1)+ ...
	     [1;-1]*y_image(1)*((d(2)-x_image(1))/(stops(1)-x_image(1))-1)));
  end
  ray_opt(:,:,2)=[1,0;d(2),1]*[1,-1/Fr(1);0,1]*ray_opt(:,:,1);
  if n(1)==-1&~reflinv,ray_opt(:,:,2)=-ray_opt(:,:,2);end
end
if ~isnan(theta0),
  y_probes(imin:imax)=w0(1)*abs(1-(x(imin:imax)-dist(1))/Fr(1))+ ...
	      theta0(1)*abs(x(imin:imax)-dist(1)*(x(imin:imax)-dist(1))/Fr(1));
else
  y_probes(imin:imax)=max(abs(repmat(apertures(1)* ...
	       (x(imin:imax)-dist(1)-x_image(1))/ ...
				   (stops(1)-x_image(1)),2,1)+ ...
		  [1;-1]*y_image(1)*((x(imin:imax)-dist(1)-x_image(1))/ ...
				     (stops(1)-x_image(1))-1)));
end
imin=imax+1;
for i=2:length(d),
  x_image(i)=1/(1/Fr(i)-1/(d(i)-x_image(i-1)));
  x_focus(i)=1/(1/Fr(i)-1/(d(i)-x_focus(i-1)));
  if ~isnan(theta0),
    theta(:,i)=theta(:,i-1)*(1-(d(i)-x_image(i-1))/Fr(i))-y_image(i-1)/Fr(i);
    if n(i)==-1&~reflinv,
      theta(:,i)=-theta(:,i);
    end
  else
    stops(i)=1/(1/Fr(i)-1/(d(i)-stops(i-1)));
    apertures(i)=abs(apertures(i-1)*Fr(i)/(d(i)-stops(i-1)-Fr(i)));
  end
  y_image(i)=-y_image(i-1)*Fr(i)/(d(i)-x_image(i-1)-Fr(i));
  if n(i)==-1&~reflinv,
    y_image(i)=-y_image(i);
  end
  y_focus(i)=abs(y_focus(i-1)*Fr(i)/(d(i)-x_focus(i-1)-Fr(i)));
  ray_im(:,1:2,i)=[1,0;x_image(i),1]*[1,-1/Fr(i);0,1]*ray_opt(:,1:2,i);
  ray_foc(:,1:2,i)=[1,0;x_focus(i),1]*[1,-1/Fr(i);0,1]*ray_opt(:,3:4,i);
  if n(i)==-1&~reflinv,
    ray_im(:,1:2,i)=-ray_im(:,1:2,i);
    ray_foc(:,1:2,i)=-ray_foc(:,1:2,i);
  end
  if i==length(d),
    imax=length(x);
  else
    imax=max(find(x<=dist(i+1)));
    if isempty(imax),imax=0;end
    if ~isnan(theta0),
      h_optics(i+1)=max(abs(y_image(i-1)*(1-d(i+1)/Fr(i))+ ...
     		theta(:,i-1)*(d(i+1)+d(i)-x_image(i-1)-d(i+1)* ...
			(d(i)-x_image(i-1))/Fr(i))));
    else
      h_optics(i+1)=max(abs(repmat(apertures(i)*(d(i+1)-x_image(i))/ ...
				      (stops(i)-x_image(i)),2,1)+ ...
		[1;-1]*y_image(i)*((d(i+1)-x_image(i))/ ...
				   (stops(i)-x_image(i))-1)));
    end
    ray_opt(:,:,i+1)=[1,0;d(i+1),1]*[1,-1/Fr(i);0,1]*ray_opt(:,:,i);
    if n(i)==-1&~reflinv,ray_opt(:,:,i+1)=-ray_opt(:,:,i+1);end
  end
  if ~isnan(theta0),
    y_probes(imin:imax)=max(abs([1;1]*y_image(i-1)* ...
				(1-(x(imin:imax)-dist(i))/Fr(i))+ ...
	       theta(:,i-1)*(x(imin:imax)-dist(i-1)-x_image(i-1)- ...
		    (d(i)-x_image(i-1))*(x(imin:imax)-dist(i))/Fr(i))));
  else
    y_probes(imin:imax)=max(abs(repmat(apertures(i)* ...
			      x(imin:imax)-dist(i)-x_image(i),2,1)/ ...
				   (stops(i)-x_image(i))+ ...
		  [1;-1]*y_image(i)*((x(imin:imax)-dist(i)-x_image(i))/ ...
				     (stops(i)-x_image(i))-1)));
  end
  imin=imax+1;
end
M_image=[y_image(1)/w0(1),(y_image(2:end)+eps)./(y_image(1:end-1)+eps)];
Mtot=prod(M_image);
y_image=abs(y_image);
M_focus=-[0,(x_focus(2:end)+eps)./(d(2:end)-x_focus(1:end-1)+eps)];
if ~reflinv,
  M_focus(n==-1)=-M_focus(n==-1);
end

% aberrations from each optic at the image created by that optic
[aberr_o,blur,OPD] = seidel_third_aberr(n,[C1;C2],[A2_1;A2_2],[A4_1;A4_2], ...
		shiftdim(ray_opt(1,1,:),1),shiftdim(ray_opt(2,1,:),1), ...
		shiftdim(ray_opt(1,2,:),1),shiftdim(ray_opt(2,2,:),1));
% aberrations from each optic at the focus created by that optic
[aberr_i,blur,OPD] = seidel_third_aberr(n,[C1;C2],[A2_1;A2_2],[A4_1;A4_2], ...
		shiftdim(ray_opt(1,3,:),1),shiftdim(ray_opt(2,3,:),1), ...
		shiftdim(ray_opt(1,4,:),1),shiftdim(ray_opt(2,4,:),1));
l_a=size(aberr_o,1);
for i=1:size(aberr_o,2),
  cum_aberr_o(:,i)=sum(aberr_o(:,1:i).* ...
		     repmat(fliplr(cumprod([1,M_image(i:-1:2)])),l_a,1),2);
end
for i=1:size(aberr_i,2),
  cum_aberr_i(:,i)=sum(aberr_i(:,1:i).* ...
		     repmat(fliplr(cumprod([1,M_focus(i:-1:2)])),l_a,1),2);
end
% cumulative OPDs
wweights=[1/4,1,0,1/2,1/2,1];
optweights=[1/4,2/3,0,1,1,1];
cum_aberr_o(end+1:end+3,:)=([2*optweights;wweights;wweights.*optweights]* ...
		abs(cum_aberr_o)).*[ones(1,size(aberr_o,2)); ...
			repmat(shiftdim(abs(ray_im(1,1,:)),1),2,1)];
cum_aberr_i(end+1:end+3,:)=([2*optweights;wweights;wweights.*optweights]* ...
		abs(cum_aberr_i)).*[ones(1,size(aberr_i,2)); ...
			repmat(shiftdim(abs(ray_foc(1,1,:)),1),2,1)];

aberr_fields={'TSC','SCC','TCC','TAC','TPC','DC','blur','OPDw','OPDo'};
aberr_o=cell2struct(num2cell(cum_aberr_o),aberr_fields,1)';
aberr_i=cell2struct(num2cell(cum_aberr_i),aberr_fields,1)';

% Complete raytracing by extending rays to image and a bit beyond
% intersections with plane located approximately one longitudinal aberration
%   beyond image plane
[l,s_nor] = ray_aspheric_intersect( ...
	     rays(:,:,end)-dist(end)-x_image(end) ...
		-2*abs(aberr_o(end).blur)/abs(ray_im(1,1,end)), ...
	     raydir(:,:,end),1,0,0,0,d(end));
raytracing(:,:,end+1)=rays(:,:,end)+raydir(:,:,end).*repmat(l,3,1);
%
raytracing=permute(raytracing,[3 1 2]);

if plotflag,
% Plots
Rarr=linspace(-1,1,100)';
fanarr=linspace(-1,1,7)';
zarr=linspace(-1,1,7)';

% prepare arrays
opt_axis=[[opt_axis';dist(end)+x_image(end)],...
	  zeros(length(opt_axis)+1,1)];
ray_a=[opt_axis([1,2:2:end],1),[0;squeeze(ray_opt(2,1,:));ray_im(2,1,end)]];
ray_p=[opt_axis([1,2:2:end],1), ...
      [w0(1);squeeze(ray_opt(2,2,:));ray_im(2,2,end)]];
lray=size(opt_axis,1)/2+1;
for i=1:length(d),
  ytemp=y_adjust(i)+h_optics(i)*Rarr-oad(i);
  if shape(i)>0&shape(i)<6,
    if abs(beta_opt(1,i))<5*eps,
      tempdraw1=alpha_opt(1,i)*ytemp.^2;
    else
      discrim=1-4*alpha_opt(1,i)*beta_opt(1,i)*ytemp.^2;
      if discrim<-5*eps,
	disp(['Error in drawing optic ',int2str(i)]);
	return
      elseif discrim<0,
        discrim=0;
      end
      tempdraw1=0.5/beta_opt(1,i)*(1-sqrt(discrim));
    end
    if n(i)~=-1,
      if abs(beta_opt(2,i))<5*eps,
        tempdraw2=alpha_opt(2,i)*ytemp.^2;
      else
        discrim=1-4*alpha_opt(2,i)*beta_opt(2,i)*ytemp.^2;
        if discrim<-5*eps,
	  disp(['Error in drawing optic ',int2str(i)]);
	  return
        elseif discrim<0,
          discrim=0;
        end
        tempdraw2=0.5/beta_opt(2,i)*(1-sqrt(discrim));
      end
    else
      tempdraw2=2*xshift(i)-tempdraw1;
    end
  else
    discrim1=1-(C1(i)*ytemp).^2;
    discrim2=1-(C2(i)*ytemp).^2;
    if discrim1<-5*eps,
      disp(['Error in drawing optic ',int2str(i)]);
      return
    elseif discrim1<0,
      discrim1=0;
    end
    if discrim2<-5*eps,
      disp(['Error in drawing optic ',int2str(i)]);
      return
    elseif discrim2<0,
      discrim2=0;
    end
    tempdraw1=C1(i)*ytemp.^2./(1+sqrt(discrim1))+ ...
              A2_1(i)*ytemp.^2+A4_1(i)*ytemp.^4;
    tempdraw2=C2(i)*ytemp.^2./(1+sqrt(discrim2))+ ...
              A2_2(i)*ytemp.^2+A4_2(i)*ytemp.^4;
  end
  opt_draw(:,:,i,1)=[tempdraw1,y_adjust(i)+h_optics(i)*Rarr]* ...
		     [cos(oaa(i)),sin(oaa(i));-sin(oaa(i)),cos(oaa(i))]+ ...
             repmat([dist(i)-(pri1(i)-xshift(i))/cos(oaa(i))-xshift(i),0], ...
		    length(Rarr),1);
  opt_draw(:,:,i,2)=[tempdraw2+thickness(i),y_adjust(i)+h_optics(i)*Rarr]* ...
		     [cos(oaa(i)),sin(oaa(i));-sin(oaa(i)),cos(oaa(i))]+ ...
             repmat([dist(i)-(pri1(i)-xshift(i))/cos(oaa(i))-xshift(i),0], ...
		    length(Rarr),1);
end
if ~isnan(theta0),
  for i=1:length(zarr),
    clear parmat
    parmat(:,:,1)=[theta0(1)*fanarr';repmat(w0(1)*zarr(i),1,length(fanarr))];
    for j=1:length(d),
      parmat(:,:,j+1)=[1,-1/Fr(j);0,1]*[1,0;d(j),1]*parmat(:,:,j);
      if n(j)==-1&~reflinv,parmat(:,:,j+1)=-parmat(:,:,j+1);end
    end
    parmat(:,:,end+1)=[1,0;x_image(end),1]*parmat(:,:,end);
    bundle(:,:,:,i)=[repmat(opt_axis([1,2:2:end],1),[1,1,length(fanarr)]), ...
		     permute(parmat(2,:,:),[3 1 2])];
  end
else
  for i=1:length(zarr),
    clear parmat
    parmat(:,:,1)=[(-aperture(1)-w0(1)*zarr(i))/stop(1)+(fanarr'-fanarr(1))* ...
		     2*aperture(1)/stop(1)/(fanarr(end)-fanarr(1)); ...
		   repmat(w0(1)*zarr(i),1,length(fanarr))];
    for j=1:length(d),
      parmat(:,:,j+1)=[1,-1/Fr(j);0,1]*[1,0;d(j),1]*parmat(:,:,j);
      if n(j)==-1&~reflinv,parmat(:,:,j+1)=-parmat(:,:,j+1);end
    end
    parmat(:,:,end+1)=[1,0;x_image(end),1]*parmat(:,:,end);
    bundle(:,:,:,i)=[repmat(opt_axis([1,2:2:end],1),[1,1,length(fanarr)]), ...
		     permute(parmat(2,:,:),[3 1 2])];
% intercepts with stop
    y_bundle_stop(i,:)=[stop(1),1]*parmat(:,:,1);
  end
end
% Rotate optical axis through off-axis optics
for i=length(d):-1:1,
  rotan=-pivot_an(i);
  rotmat=[cos(rotan),-sin(rotan);sin(rotan),cos(rotan)];
  ray_a(i+2:end,:)=(ray_a(i+2:end,:)+ ...
		    repmat(origin(1:2,i)',size(ray_a,1)-i-1,1))*rotmat;
  ray_p(i+2:end,:)=(ray_p(i+2:end,:)+ ...
		    repmat(origin(1:2,i)',size(ray_p,1)-i-1,1))*rotmat;
  for k=1:length(zarr),
   for j=1:length(fanarr),
    bundle(i+1:end,:,j,k)=(bundle(i+1:end,:,j,k)+ ...
		   repmat(origin(1:2,i)',size(bundle,1)-i,1))*rotmat;
   end
  end
  for k=1:size(raytracing,3),
    raytracing(2*i+1:end,1:2,k)=(raytracing(2*i+1:end,1:2,k)+ ...
		   repmat(origin(1:2,i)',size(raytracing,1)-2*i,1))*rotmat;
  end
  for k=1:2,
   for j=i+1:length(d),
    opt_draw(:,:,j,k)=(opt_draw(:,:,j,k)+ ...
		   repmat(origin(1:2,i)',size(opt_draw,1),1))*rotmat;
   end
  end
  opt_axis(2*i+1:end,:)=(opt_axis(2*i+1:end,:)+ ...
		   repmat(origin(1:2,i)',size(opt_axis,1)-2*i,1))*rotmat;
end
% Axial and principal rays
fc1=figure(19631);
fpos=get(fc1,'Position');
set (fc1,'Position',[10 20 fpos(3:4)])
% object
plot([0;0],[0;w0(1)],'LineWidth',6);
hold on
% axial ray
plot(ray_a(:,1),ray_a(:,2),'k-');
% principal ray
plot(ray_p(:,1),ray_p(:,2),'k-');
% image
plot([ray_a(end,1);ray_p(end,1)],[ray_a(end,2);ray_p(end,2)],'LineWidth',6);
% optics
for i=1:length(d),
  plot(opt_draw(:,1,i,1),opt_draw(:,2,i,1),'m');
  if n(i)~=-1,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m');
  elseif reflinv,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m--');
  end
end
if isnan(theta0),
% stop
  plot(stop(1)*[1;1],aperture(1)*[1;1.5],'k','LineWidth',8);
  plot(stop(1)*[1;1],-aperture(1)*[1;1.5],'k','LineWidth',8);
% rays to stop if behind object or if past first optic
  if stop(1)<0,
    plot([0;stop(1)],[0;aperture(1)],'k--');
    plot([0;stop(1)],[w0(1);0],'k--');
  end
  if stop(1)>d(1),
    plot([ray_a(2,1);stop(1)],[ray_a(2,2);aperture(1)],'k--');
    plot([ray_p(2,1);stop(1)],[ray_p(2,2);0],'k--');
  end
end
% optical axis
for i=1:size(opt_axis,1)/2,
  plot(opt_axis(2*i-1:2*i,1),opt_axis(2*i-1:2*i,2),':');
end
if isnan(theta0)&stop(1)<0,
  plot([opt_axis(1,1);stop(1)],opt_axis(1,2)*[1;-1],':');
end

title('Axial and principal rays')
xlabel('x (m)')
ylabel('y (m)')

% Paraxial ray bundle
fc2=figure(19632);
fpos=get(fc2,'Position');
set (fc2,'Position',[620 530 fpos(3:4)])
% object
plot([0;0],[-w0(1);w0(1)],'LineWidth',6)
hold on
% ray bundle
for i=1:length(zarr),
  plot(squeeze(bundle(:,1,:,i)),squeeze(bundle(:,2,:,i)));
end
if isnan(theta0),
% rays to stop if behind object or if past first optic
  if stop(1)<0,
    for i=1:length(zarr),
      plot([0;stop(1)],[squeeze(bundle(1,2,:,i))';y_bundle_stop(i,:)],'--');
    end
  end
  if stop(1)>d(1),
    for i=1:length(zarr),
      plot([opt_axis(2,1);stop(1)], ...
	   [squeeze(bundle(2,2,:,i))';y_bundle_stop(i,:)],'--');
    end
  end
% stop
  plot(stop(1)*[1;1],aperture(1)*[1;1.5],'k','LineWidth',8);
  plot(stop(1)*[1;1],-aperture(1)*[1;1.5],'k','LineWidth',8);
end
% image
plot([ray_p(end,1);2*ray_a(end,1)-ray_p(end,1)], ...
     [ray_p(end,2);2*ray_a(end,2)-ray_p(end,2)],'LineWidth',6);
% optics
for i=1:length(d),
  plot(opt_draw(:,1,i,1),opt_draw(:,2,i,1),'m');
  if n(i)~=-1,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m');
  elseif reflinv,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m--');
  end
end
% optical axis
for i=1:size(opt_axis,1)/2,
  plot(opt_axis(2*i-1:2*i,1),opt_axis(2*i-1:2*i,2),':');
end

title('Paraxial ray bundle')
xlabel('x (m)')
ylabel('y (m)')

% Detailed ray tracing
fc3=figure(19633);
fpos=get(fc3,'Position');
set (fc3,'Position',[650 10 fpos(3:4)])
% object
plot([0;0],[-w0(1);w0(1)],'LineWidth',6)
hold on
% rays
t1=find(n==-1);
j=1;
for i=t1,
  plot(squeeze(raytracing(j:2*i,1,1:nray3)), ...
       squeeze(raytracing(j:2*i,2,1:nray3)));
  j=2*i+1;
end
plot(squeeze(raytracing(j:end,1,1:nray3)),squeeze(raytracing(j:end,2,1:nray3)));
if isnan(theta0),
% rays to stop if behind object or if past first optic
  if stop(1)<0,
    plot([0;stop(1)],[squeeze(raytracing(1,2,1:nray3))'; ...
		              y_ray_stop(1:nray3)],'--');
  end
  if stop(1)>d(1),
    plot([opt_axis(2,1);stop(1)], ...
	 [squeeze(raytracing(2,2,1:nray3))';y_ray_stop(1:nray3)],'--');
  end
% stop
  plot(stop(1)*[1;1],aperture(1)*[1;1.5],'k','LineWidth',8);
  plot(stop(1)*[1;1],-aperture(1)*[1;1.5],'k','LineWidth',8);
end
% image
plot([ray_p(end,1);2*ray_a(end,1)-ray_p(end,1)], ...
     [ray_p(end,2);2*ray_a(end,2)-ray_p(end,2)],'LineWidth',6);
% optics
for i=1:length(d),
  plot(opt_draw(:,1,i,1),opt_draw(:,2,i,1),'m');
  if n(i)~=-1,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m');
  elseif reflinv,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m--');
  end
end
% optical axis
for i=1:size(opt_axis,1)/2,
  plot(opt_axis(2*i-1:2*i,1),opt_axis(2*i-1:2*i,2),':');
end

title('Ray tracing')
xlabel('x (m)')
ylabel('y (m)')

% Detailed ray tracing (single bundle)
fc4=figure(19634);
fpos=get(fc4,'Position');
set (fc4,'Position',[10 520 fpos(3:4)])
% object
plot([0;0],[-w0(1);w0(1)],'LineWidth',6)
hold on
% rays
t1=find(n==-1);
j=1;
for i=t1,
  plot(squeeze(raytracing(j:2*i,1,nray3+1:end)), ...
       squeeze(raytracing(j:2*i,2,nray3+1:end)));
  j=2*i+1;
end
plot(squeeze(raytracing(j:end,1,nray3+1:end)), ...
     squeeze(raytracing(j:end,2,nray3+1:end)));
if isnan(theta0),
% rays to stop if behind object or if past first optic
  if stop(1)<0,
    plot([0;stop(1)],[squeeze(raytracing(1,2,nray3+1:end))'; ...
			      y_ray_stop(nray3+1:end)],'--');
  end
  if stop(1)>d(1),
    plot([opt_axis(2,1);stop(1)], ...
	 [squeeze(raytracing(2,2,nray3+1:end))';y_ray_stop(nray3+1:end)],'--');
  end
% stop
  plot(stop(1)*[1;1],aperture(1)*[1;1.5],'k','LineWidth',8);
  plot(stop(1)*[1;1],-aperture(1)*[1;1.5],'k','LineWidth',8);
end
% image
plot([ray_p(end,1);2*ray_a(end,1)-ray_p(end,1)], ...
     [ray_p(end,2);2*ray_a(end,2)-ray_p(end,2)],'LineWidth',6);
% optics
for i=1:length(d),
  plot(opt_draw(:,1,i,1),opt_draw(:,2,i,1),'m');
  if n(i)~=-1,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m');
  elseif reflinv,
    plot(opt_draw(:,1,i,2),opt_draw(:,2,i,2),'m--');
  end
end
% optical axis
for i=1:size(opt_axis,1)/2,
  plot(opt_axis(2*i-1:2*i,1),opt_axis(2*i-1:2*i,2),':');
end
% limit plot to a reasonable size to the left
xlim=get(gca,'xlim');
xlim(1)=max(min(stop(1),-2*abs(xlim(2))),xlim(1));
set(gca,'xlim',xlim);
if mray4>1,
  title(['Ray tracing for parallel rays (', ...
		num2str(sqrt(1-raydir_y^2-raydir_z^2)),',', ...
		num2str(raydir_y),',',num2str(raydir_z),')'])
else
  title(['Ray tracing for point x=',num2str(raysource_x*w0(1)), ...
			   ', y=',num2str(raysource_y*w0(1)), ...
			   ', z=',num2str(raysource_z*w0(1))])
end
xlabel('x (m)')
ylabel('y (m)')
end

% reorder probes points
x=x(probesorder);
y_probes=y_probes(probesorder);

% update input parameters
for i=1:l_in,
 for k=1:length(input_data{i}),
  for fi=input_fields{i},
   input_data{i}=setfield(input_data{i},{k},fi{1},eval([fi{1},'(k)']));
  end
 end
end

% fold outputs into output parameters
output_variables={'image','focus','optics','probes','Mtot'};

output_fields{1}={'x','y','u','up','M','aberrations'};
output_fields{2}={'x','y','M','aberrations'};
output_fields{3}={'h','u','y','up','yp'};
output_fields{4}={'y'};
output_fields{5}={[]};

if ~isempty(x_image),
 image=cell2struct([num2cell([x_image;y_image; ...
		    reshape(ray_im(1,1:2,:),2,length(optics));M_image]); ...
		    num2cell(aberr_o)], ...
		   output_fields{strcmp(output_variables,'image')},1);
end
if ~isempty(x_focus),
 focus=cell2struct([num2cell([x_focus;y_focus;M_focus]);num2cell(aberr_i)], ...
		  output_fields{strcmp(output_variables,'focus')},1);
end
if ~isempty(h_optics),
 is=find(strcmp(input_variables,'optics'));
 cell_data{is}=struct2cell(orderfields(input_data{is}));
 optics=cell2struct([cell_data{is}; ...
          num2cell([h_optics;reshape(ray_opt(:,1:2,:),4,length(optics))])], ...
	[input_fields{is},output_fields{strcmp(output_variables,'optics')}],1);
end
if ~isempty(y_probes),
 is=find(strcmp(input_variables,'probes'));
 cell_data{is}=struct2cell(orderfields(input_data{is}));
 probes=cell2struct([cell_data{is};num2cell(y_probes)], ...
	[input_fields{is},output_fields{strcmp(output_variables,'probes')}],1);
end

%
%

function [aberrations,blur,OPD] = seidel_third_aberr(n,C,A2,A4,u,y,up,yp)

% all calculations for each optic are done at image generated *by that optic*

% n		= indexes of refraction (row vector)
% C (m-1)	= curvatures (2 rows)
% A2 (m-1)	= aspheric coefficients (2 rows)
% A4 (m-3)	= aspheric coefficients (2 rows)
% u (rad) 	= angles of axial (marginal) rays (row vector)
%		    convention: positive when shifted counterclockwise from axis
% y (m) 	= heights of axial (marginal) rays (row vector)
% up (rad) 	= angles of principal (chief) rays (row vector)
%		    convention: positive when shifted counterclockwise from axis
% yp (m) 	= heights of principal (chief) rays (row vector)

% effective curvature
Ce=C+2*A2;
Ctot=-diff(Ce);

% power
phi=Ctot.*(n-1);

v=u./y;
Q=yp./y;

iu=find(u);inu=find(~u);
% object at finite distance
p(iu)=y(iu)./u(iu);
q(iu)=1./(phi(iu)-1./p(iu));
% object at infinity (focusing optic)
p(inu)=Inf;
q(inu)=1./phi(inu);

uprime=u-y.*phi;
upprime=up-yp.*phi;
h_im=q.*upprime+yp;

% these are the equations for a thin lens, written here for reference; here,
%   however, we calculate refraction from individual surfaces to cover reflectors
%   and aspherics
%G1=n.^2.*(n-1)/2;
%G2=(2*n+1).*(n-1)/2;
%G3=(3*n+1).*(n-1)/2;
%G4=(n+2).*(n-1)/2./n;
%G5=2*(n+1).*(n-1)./n;
%G6=(3*n+2).*(n-1)/2./n;
%G7=(2*n+1).*(n-1)/2./n;
%G8=n.*(n-1)/2;
%
%TSC=y.^4.*(G1.*Ctot.^3-G2.*Ctot.^2.*Ce(1,:)-G3.*Ctot.^2.*v ...
%	    +G4.*Ctot.*Ce(1,:).^2+G5.*Ctot.*Ce(1,:).*v+G6.*Ctot.*v.^2)./uprime;
%SCC=-h_im.*y.^2.*(G5.*Ctot.*Ce(1,:)/4+G7.*Ctot.*v-G8.*Ctot.^2);
%TAC=h_im.^2.*phi.*uprime/2;
%TPC=TAC./n;
%DC=0;
%
%DC=DC+Q.*(TPC+3*TAC)+3*Q.^2.*SCC+Q.^3.*TSC;
%TAC=TAC+2*Q.*SCC+Q.^2.*TSC;
%SCC=SCC+Q.*TSC;
%TCC=3*SCC;
%

% optical invariant
inv=yp.*u-y.*up;inv=inv(1);
% angles of incidence on first surface
i=y.*C(1,:)+u;
ip=yp.*C(1,:)+up;
% u angles after exiting both surfaces
u_after=[u-i+i./n;uprime];
up_after=[up-ip+ip./n;upprime];
% angles of incidence on both surfaces
i=[i;y.*C(2,:)+u_after(1,:)];
ip=[ip;yp.*C(2,:)+up_after(1,:)];
% optical parameters
B=[1-1./n;n.*(1-n)].*repmat(y,2,1).*(u_after+i)/2/inv;
Bp=[1-1./n;n.*(1-n)].*repmat(yp,2,1).*(up_after+ip)/2/inv;
K=A4-A2.*(4*A2.^2+6*C.*A2+3*C.^2)/4;
W=4*K.*[n-1;1-n].*repmat(h_im,2,1)/inv;

% TSC, SCC, TCC, TAC, TPC, DC
aberrations(1,:)=sum(B.*i.^2.*repmat(h_im,2,1)+W.*repmat(y.^4,2,1));
aberrations(2,:)=sum(B.*i.*ip.*repmat(h_im,2,1)+W.*repmat(y.^3.*yp,2,1));
aberrations(3,:)=3*aberrations(2,:);
aberrations(4,:)=sum(B.*ip.^2.*repmat(h_im,2,1)+W.*repmat(y.^2.*yp.^2,2,1));
aberrations(5,:)=inv/2*sum([n-1;1-n].*Ce.*repmat(h_im./n,2,1));
aberrations(6,:)=sum((Bp.*i.*ip+(up_after.^2-[up;up_after(1,:)].^2)/2).* ...
		 repmat(h_im,2,1)+W.*repmat(y.*yp.^3,2,1));
% optimized blur
blur=[1/4,2/3,0,1,1,1]*abs(aberrations);
% worst-case and optimized OPD
OPD(1,:)=([1/4,1,0,1/2,1/2,1]*abs(aberrations)).*abs(uprime);
OPD(2,:)=([1/16,2/3,0,1/2,1/2,1]*abs(aberrations)).*abs(uprime);


function [l,s_nor] = ray_aspheric_intersect(R0,r,shape,C,A2,A4,l0)

% intersections of a set of rays (origin R0, direction vector r) with an 
%   aspheric surface, either specified as an exact conic or with a 
%   near-spherical approximation (shape=0)
%	conic:	shape=1 or 2: spherical (curvature C)
%		shape=3, 4 or 5: paraboloidal, ellipsoidal or hyperboloidal
%		  (equation x=alpha*s^2+beta*x^2, C=2*alpha, 
%		  A4=-alpha^3+alpha^2*beta)
%       coordinate system: x along axis of optic, origin at vertex
%	output:
%		l = distance along ray of intersection point
%		s_nor = direction vector of normal to the surface

if shape>0&shape<6,
% conic
  alpha_opt=C/2;
  if alpha_opt,
    beta_opt=A4/alpha_opt^2+alpha_opt;
  else
    beta_opt=0;
  end
  for i=1:size(r,2),
    l(i)=fzero(@(l) conic_intersect(l,R0(:,i),r(:,i),alpha_opt,beta_opt),l0);
  end
  x=R0+repmat(l,3,1).*r;
  s=sqrt(x(2,:).^2+x(3,:).^2);
  s(~s)=eps;  
%  s_nor=[ones(1,length(s)); ...
%	 repmat((-2*sqrt(discrim)./(1-2*beta_opt*x(1,:)))./s,2,1).*x(2:3,:)];
  s_nor=[ones(1,length(s)); ...
	 repmat(-2*alpha_opt./(1-2*beta_opt*x(1,:)),2,1).*x(2:3,:)];
  s_nor=s_nor./repmat(sqrt(sum(s_nor.^2)),3,1);
else
% spherical or near-spherical surface
  for i=1:size(r,2),
    l(i)=fzero(@(l) near_spher_intersect(l,R0(:,i),r(:,i),C,A2,A4),l0);
  end
  x=R0+repmat(l,3,1).*r;
  s=sqrt(x(2,:).^2+x(3,:).^2);
  s(~s)=eps;  
  discrim=1-C^2*s.^2;
  if any(discrim<-5*eps),
    disp('Intersection with optic not found');
    return
  elseif any(discrim<0),
    discrim(discrim<0)=0;
  end
  s_nor=[ones(1,length(s)); ...
	 repmat((-(C*s./sqrt(discrim)+2*A2*s+4*A4*s.^3))./s,2,1).*x(2:3,:)];
  s_nor=s_nor./repmat(sqrt(sum(s_nor.^2)),3,1);
end

function d = conic_intersect(l,R0,r,alpha_opt,beta_opt)

s=sqrt(sum((R0(2:3)+l*r(2:3)).^2));
d=R0(1)+l*r(1)-alpha_opt*s^2-beta_opt*(R0(1)+l*r(1))^2;

function d = near_spher_intersect(l,R0,r,C,A2,A4)

x=R0(1)+l*r(1);
s=sqrt(sum((R0(2:3)+l*r(2:3)).^2));
if abs(C*s)>=1,
  d=x-C*s^2-A2*s^2-A4*s^4;
else
  discrim=1-(C*s)^2;
  if discrim<-5*eps,
    disp('Intersection with optic not found');
    return
  elseif discrim<0,
    discrim=0;
  end
  d=x-C*s^2/(1+sqrt(discrim))-A2*s^2-A4*s^4;
end

function d = aspheric_skew_intersects(x,yp,x0,y0,theta,C,A2,A4)

s=-(x-x0)*tan(theta)+yp/cos(theta)-y0;
if abs(C*s)>=1,
  d=x-C*s^2-A2*s^2-A4*s^4;
else
  discrim=1-(C*s)^2;
  if discrim<-5*eps,
    disp('Intersection with optic not found');
    return
  elseif discrim<0,
    discrim=0;
  end
  d=x-C*s^2/(1+sqrt(discrim))-A2*s^2-A4*s^4;
end
