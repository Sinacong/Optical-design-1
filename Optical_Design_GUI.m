function optical_design_GUI(varargin)

%Matlab graphical interface to file /mac/coda/matlab/optical_design.m
%Just call the routine without any input for normal usage.
%Default parameters are also indicated by [].
%Call the routine with the second argument defined as the ratio of the resulting window divided by the ratio of the screen. This number must be between 0 and 1.
%
%AM 12/11/07

global A
what='begin';
def=0.5;
a=def;
if nargin==1
	if ~isempty(varargin{1});
		what=varargin{1};
	end
end
if nargin==2
	if ~isempty(varargin{1});
		what=varargin{1};
	end
	if ~isempty(varargin{2});
		a=varargin{2};
	end
end
if ~isstr(what)
	disp('The first argument must be a string or empty matrix, please.')
	return
end
if ~isnumeric(a)
	disp('The window size must be a number, taking default value')
	a=def;
end
if or(a<=0,a>1)
	disp('Window size must be larger than 0 and smaller or equal to 1, taking default value')
	a=def;
end

switch what
	
	case 'begin'
	
		A=struct('beam',[],'optics',[],'probes',[]);
		set(0,'DefaultUicontrolUnits','normalized')
		scr=get(0,'screensize');
		f(1)=(1-a)/2;
		f(2)=(1-a)/2;
		f(3)=a;
		f(4)=a/2;
		h=figure('units','normalized','position',[f(1) f(2) f(3) f(4)]);
		set(h,'name','Matlab interface: Optical design calculations','numbertitle','off');

		%Beam parameters
		uicontrol('style','text','position',[0.02 0.90 0.4 0.1],'string','Beam parameters');
		uicontrol('style','text','position',[0.02 0.8 0.3 0.1],'string','Half width [m]');
		uicontrol('style','text','position',[0.02 0.7 0.3 0.1],'string','Half divergence [rad]');
		uicontrol('style','text','position',[0.02 0.6 0.3 0.1],'string','Aperture stop [m]');
		uicontrol('style','text','position',[0.02 0.5 0.3 0.1],'string','Obj-stop distance [m]');

		uicontrol('style','edit','position',[0.3 0.8 0.1 0.1],'tag','bhw','string','0.035');
		uicontrol('style','edit','position',[0.3 0.7 0.1 0.1],'tag','bhdiv','string','0.01');
		uicontrol('style','edit','position',[0.3 0.6 0.1 0.1],'tag','ahs','string','0.045');
		uicontrol('style','edit','position',[0.3 0.5 0.1 0.1],'tag','objstopdist','string','1');
		uicontrol('style','popupmenu','position',[0.02 0.4 0.4 0.1],'tag','what','string','Half divergence|Aperture stop');


		%Optics parameters
		number=mat2cell(1:10,1,10);
		uicontrol('style','text','position',[0.5 0.9 0.2 0.1],'string','Optics parameters');

		uicontrol('style','text','position',[0.5 0.8 0.2 0.1]','string','Number of optics');
		uicontrol('position',[0.5 0.6 0.1 0.1],'callback','optical_design_GUI optics');
		uicontrol('style','listbox','position',[0.6 0.6 0.1 0.1],'string',number,'tag','list');
		
		%Probes
		uicontrol('style','text','position',[0.75 0.9 0.2 0.1],'string','Probes');

		uicontrol('style','text','position',[0.75 0.8 0.2 0.1]','string','Number of probes');
		uicontrol('style','edit','position',[0.75 0.7 0.1 0.1],'string','2','tag','numpr');
		uicontrol('position',[0.75 0.6 0.1 0.1],'callback','optical_design_GUI probes');
		
		%Calculate
		uicontrol('position',[0.8 0.2 0.1 0.1],'string','calculate','callback','optical_design_GUI calculate')
		
	case 'optics'

		n=get(findobj(gcbf,'style','listbox','tag','list'),'value');
		for i=1:n
			hhh=figure;
			set(0,'DefaultUicontrolUnits','normalized');
			set(hhh,'name','Optical design calculations, optics parameters','numbertitle','off');
			uicontrol('position',[0.1 0.94 0.3 0.05],'style','text','string','Focal length [m]');
			uicontrol('position',[0.1 0.85 0.3 0.05],'style','text','string','Distance from the previous one [m]');
			uicontrol('position',[0.1 0.75 0.3 0.05],'style','text','string','Curvature of entry surface [m-1]');
			uicontrol('position',[0.1 0.65 0.3 0.05],'style','text','string','Curvature of exit surface [m-1]');
			uicontrol('position',[0.1 0.55 0.3 0.05],'style','text','string','Index of refraction');
			uicontrol('position',[0.1 0.45 0.3 0.05],'style','text','string','A2_1 [m-1]');	
			uicontrol('position',[0.1 0.35 0.3 0.05],'style','text','string','A2_2 [m-1]');
			uicontrol('position',[0.1 0.25 0.3 0.05],'style','text','string','A4_1 [m-1]');		
			uicontrol('position',[0.1 0.15 0.3 0.05],'style','text','string','A4_2 [m-3]');			
			uicontrol('position',[0.1 0.05 0.3 0.05],'style','text','string','Shape');	
			uicontrol('position',[0.75 0.5 0.2 0.05],'style','text','string',strcat(['Element',' ',num2str(i)]));
			
			uicontrol('position',[0.4 0.94 0.3 0.05],'style','edit','string','1','tag','focal');
			uicontrol('position',[0.4 0.85 0.3 0.05],'style','edit','string','1','tag','dist');
			uicontrol('position',[0.4 0.75 0.3 0.05],'style','edit','string','0','tag','C1');
			uicontrol('position',[0.4 0.65 0.3 0.05],'style','edit','string','0','tag','C2');
			uicontrol('position',[0.4 0.55 0.3 0.05],'style','edit','string','2.54','tag','refr_ind');
			uicontrol('position',[0.4 0.45 0.3 0.05],'style','edit','string','0','tag','A21');	
			uicontrol('position',[0.4 0.35 0.3 0.05],'style','edit','string','0','tag','A22');
			uicontrol('position',[0.4 0.25 0.3 0.05],'style','edit','string','0','tag','A41');		
			uicontrol('position',[0.4 0.15 0.3 0.05],'style','edit','string','0','tag','A42');			
			uicontrol('position',[0.4 0.05 0.3 0.05],'style','popupmenu','string',...
			'generic|convex-plano lens|spherical mirror|paraboloidal mirror','tag','shape');
			uicontrol('position',[0.8 0.15 0.1 0.05],'string','Ok','callback','optical_design_GUI takeoptics',...
			'userdata',[n i]);
		end					
		
	case 'takeoptics'
	
		nn=get(gcbo,'userdata');
		A.optics(nn(2)).A2_1=str2num(get(findobj(gcbf,'style','edit','tag','A21'),'string'));
		A.optics(nn(2)).A2_2=str2num(get(findobj(gcbf,'style','edit','tag','A22'),'string'));
		A.optics(nn(2)).A4_1=str2num(get(findobj(gcbf,'style','edit','tag','A41'),'string'));
		A.optics(nn(2)).A4_2=str2num(get(findobj(gcbf,'style','edit','tag','A42'),'string'));
		A.optics(nn(2)).C1=str2num(get(findobj(gcbf,'style','edit','tag','C1'),'string'));
		A.optics(nn(2)).C2=str2num(get(findobj(gcbf,'style','edit','tag','C2'),'string'));
		A.optics(nn(2)).F=str2num(get(findobj(gcbf,'style','edit','tag','focal'),'string'));							
		A.optics(nn(2)).d=str2num(get(findobj(gcbf,'style','edit','tag','dist'),'string'));
		A.optics(nn(2)).n=str2num(get(findobj(gcbf,'style','edit','tag','refr_ind'),'string'));
		A.optics(nn(2)).shape=get(findobj(gcbf,'style','popupmenu','tag','shape'),'value');
		close(gcbf)
		
	case 'probes'
	
		numpr=str2num(get(findobj(gcbf,'style','edit','tag','numpr'),'string'));
		set(0,'DefaultUicontrolUnits','normalized')
		scr=get(0,'screensize');
		hh=figure('units','normalized','position',[0.3 0.3 0.7 0.1]);
		set(hh,'name','Matlab interface: Optical design calculations, probes parameters','numbertitle','off');
		uicontrol('position',[0.05 0.55 0.8 0.4],'style','text','string','Probes positions');
		uicontrol('position',[0.05 0.05 0.8 0.4],'tag','probes','style','edit','string','');
		uicontrol('position',[0.85 0.05 0.1 0.4],'string','ok','callback',...
		'optical_design_GUI takeprobes','userdata',numpr);

	case 'takeprobes'
	
		num=get(gco,'userdata');
		probloc=str2num(get(findobj(gcbf,'style','edit','tag','probes'),'string'));
		if length(probloc)~=num
			disp('The number of probe locations does not correspond to the number of declared probes.');
			return
		end
		for i=1:num
			A.probes(i).x=probloc(i);
		end
		close(gcbf)
		
	case 'calculate'
	
		A.beam.stop=str2num(get(findobj(gcbf,'style','edit','tag','objstopdist'),'string'));
		A.beam.theta0=str2num(get(findobj(gcbf,'style','edit','tag','bhdiv'),'string'));
		A.beam.w0=str2num(get(findobj(gcbf,'style','edit','tag','bhw'),'string'));
		A.beam.aperture=str2num(get(findobj(gcbf,'style','edit','tag','ahs'),'string'));
		what=get(findobj(gcbf,'tag','what','style','popupmenu'),'value');
		if what==2
			A.beam.theta0=[];
		end
		save('latest_optical_design.mat','A');
		[B.image,B.focus,B.optics,B.probes,B.Mtot]=optical_design(A.beam,A.optics,A.probes);
		assignin('base','B',B);
end

