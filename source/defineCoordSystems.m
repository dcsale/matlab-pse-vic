function [Coord] = defineCoordSystems(BLADE, azim)

%% Global coordinate system
%  notes: this is the reference coordinate system, it does not ever move
% origin: (0,0,0), the apex of the cone of rotation
% +x dir: points in the downstream direction of the nominal free stream 
% +y dir: orthogonal to global x-z axes, forms a right handed coordinate system
% +z dir: points vertically upward, opposite the direction of gravity

RZ1 = [   cosd(BLADE.YAW),  -sind(BLADE.YAW),                0;
          sind(BLADE.YAW),   cosd(BLADE.YAW),                0;
                        0,                 0,                1];  
                          
% basis vectors of the global coord. sys.
% uG = [1; 0; 0];
% vG = [0; 1; 0];
% wG = [0; 0; 1];
uG = RZ1(:,1);
vG = RZ1(:,2);
wG = RZ1(:,3);

% origin
% OG = [0; 0; 0];
OG = [0; 0; BLADE.HUB_HT];

%% Shaft coordinate system
%  notes: this system rotates with the rotor
% origin: same as the global coordinate system
% +x dir: points in the downstream direction of the possibly tilted shaft
% +y dir: orthogonal to shaft x-z axes, forms a right handed coordinate system
% +z dir: points in the direction of blade azimuth

% c1      = cosd(SHAFT_TILT);
% c2      = cosd(AZIM(n));
% s1      = sind(SHAFT_TILT);
% s2      = sind(AZIM(n));
% Tran_GS = [ c1, s1*s2, c2*s1;
%              0,    c2,   -s2;
%            -s1, c1*s2, c1*c2];

RZ1 = [   cosd(BLADE.YAW),  -sind(BLADE.YAW),                0;
          sind(BLADE.YAW),   cosd(BLADE.YAW),                0;
                        0,                 0,                1]; 
                    
RY1 = [ cosd(BLADE.SHAFT_TILT),                0,  sind(BLADE.SHAFT_TILT);
                             0,                1,                       0;
       -sind(BLADE.SHAFT_TILT),                0,  cosd(BLADE.SHAFT_TILT)];
RX2 = [                      1,                0,                       0;
                             0,       cosd(azim),             -sind(azim);
                             0,       sind(azim),             cosd(azim)];
% Rot_GS = RY1*RX2;
Rot_GS = RZ1*RY1*RX2;


% basis vectors of the shaft coord. sys. w.r.t. global coord. sys.
uS = Rot_GS(:,1);
vS = Rot_GS(:,2);
wS = Rot_GS(:,3);

% origin of the shaft coord. sys.
OS = OG;

%% Blade coordinate system
%  notes: this system rotates with the rotor, and pitches with the blade pitch
% origin: intersection of the blade pitch axis and the blade root
% +x dir: points in the direction of the trailing edge parallel w/ the chordline as if there was zero pre-twist of the blade
% +y dir: orthogonal to blade x-z axes, forms a right handed coordinate system
% +z dir: points towards the blade tip, parallel to the blade pitch axis

RZ1 = [   cosd(BLADE.YAW),  -sind(BLADE.YAW),                0;
          sind(BLADE.YAW),   cosd(BLADE.YAW),                0;
                        0,                 0,                1]; 
                    
RY1 = [  cosd(BLADE.SHAFT_TILT),                 0, sind(BLADE.SHAFT_TILT);
                              0,                 1,                      0;
        -sind(BLADE.SHAFT_TILT),                 0, cosd(BLADE.SHAFT_TILT)];
    
RX2 = [                       1,                 0,                     0;
                              0,        cosd(azim),           -sind(azim);
                              0,        sind(azim),            cosd(azim)];
                          
RY3 = [    cosd(BLADE.PRE_CONE),                 0,   sind(BLADE.PRE_CONE);
                              0,                 1,                      0;
          -sind(BLADE.PRE_CONE),                 0,   cosd(BLADE.PRE_CONE)];
      
RZ4 = [   cosd(BLADE.BLD_PITCH),  -sind(BLADE.BLD_PITCH),                0;
          sind(BLADE.BLD_PITCH),   cosd(BLADE.BLD_PITCH),                0;
                              0,                       0,                1];     
                          
% RX5 = [                       1,                 0,                     0;
%                               0,        cosd(BLADE.TEETER),           -sind(BLADE.TEETER);
%                               0,        sind(BLADE.TEETER),            cosd(BLADE.TEETER)];
% RY5 = [    cosd(BLADE.TEETER),                 0,   sind(BLADE.TEETER);
%                             0,                 1,                    0;
%           -sind(BLADE.TEETER),                 0,   cosd(BLADE.TEETER)];
      
% Rot_GB = RY1*RX2*RY3*RZ4;
Rot_GB = RZ1*RY1*RX2*RY3*RZ4;

% basis vectors of the blade coord. sys. w.r.t. global coord. sys.
uB = -Rot_GB(:,2); %Rot_GB(:,1);
vB =  Rot_GB(:,1);  %Rot_GB(:,2);
wB =  Rot_GB(:,3);

% origin of the blade coord. sys.
OB = OG + BLADE.HUB_DIA/2.*wB;
TB = OG + BLADE.ROTOR_DIA/2.*wB;

% additional shift of the origin for VAWTs (this is a HACK) -- what do I mean to hack? well, does the code match the documentation?
if BLADE.TEETER ~= 0
    OB = OB + BLADE.ROTOR_DIA/2.*vB;
    TB = TB + BLADE.ROTOR_DIA/2.*vB;
end

%% Transformation matrices
% from global to shaft coordinate system
Tran_GS = [dot(uG,uS), dot(vG,uS), dot(wG,uS); ...
           dot(uG,vS), dot(vG,vS), dot(wG,vS); ...
           dot(uG,wS), dot(vG,wS), dot(wG,wS)];       
% from shaft to global coordinate system
Tran_SG = Tran_GS';   

% from shaft to blade coordinate system
Tran_SB = [dot(uS,uB), dot(vS,uB), dot(wS,uB); ...
           dot(uS,vB), dot(vS,vB), dot(wS,vB); ...
           dot(uS,wB), dot(vS,wB), dot(wS,wB)];       
% from blade to shaft coordinate system
Tran_BS = Tran_SB';  

% from global to blade coordinate system
Tran_GB = [dot(uG,uB), dot(vG,uB), dot(wG,uB); ...
           dot(uG,vB), dot(vG,vB), dot(wG,vB); ...
           dot(uG,wB), dot(vG,wB), dot(wG,wB)];       
% from blade to global coordinate system
Tran_BG = Tran_GB';  

%% Collect the output
% the global coordinate system
Coord.Origin.Global = OG;
Coord.Tip.Global = [];
Coord.U.Global = uG;
Coord.V.Global = vG;
Coord.W.Global = wG;
% the shaft coordinate system
Coord.Origin.Shaft = OS;
Coord.Tip.Shaft = [];
Coord.U.Shaft = uS;
Coord.V.Shaft = vS;
Coord.W.Shaft = wS;
% the blade coordinate system
Coord.Origin.Blade = OB;
Coord.Tip.Blade    = TB;
Coord.U.Blade = uB;
Coord.V.Blade = vB;
Coord.W.Blade = wB;

% the transformation matrices
Coord.Tran.GS = Tran_GS;
Coord.Tran.SG = Tran_SG;

Coord.Tran.SB = Tran_SB;
Coord.Tran.BS = Tran_BS;

Coord.Tran.GB = Tran_GB;
Coord.Tran.BG = Tran_BG;

%% plot the coordinate systems to verify
% clf
% hold on
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis([-3 3 -3 3 -3 3])
% view(-60,20)
% 
% % global
% plot3(OG(1), OG(2), OG(3), '.k')
% plot3([OG(1) OG(1)+uG(1)], [OG(2) OG(2)+uG(2)], [OG(3) OG(3)+uG(3)],'-r')
% plot3([OG(1) OG(1)+vG(1)], [OG(2) OG(2)+vG(2)], [OG(3) OG(3)+vG(3)],'-b')
% plot3([OG(1) OG(1)+wG(1)], [OG(2) OG(2)+wG(2)], [OG(3) OG(3)+wG(3)],'-g')
% % shaft
% plot3(OS(1), OS(2), OS(3), '.k')
% plot3([OS(1) OS(1)+uS(1)], [OS(2) OS(2)+uS(2)], [OS(3) OS(3)+uS(3)],':r')
% plot3([OS(1) OS(1)+vS(1)], [OS(2) OS(2)+vS(2)], [OS(3) OS(3)+vS(3)],':b')
% plot3([OS(1) OS(1)+wS(1)], [OS(2) OS(2)+wS(2)], [OS(3) OS(3)+wS(3)],':g')
% % blade
% plot3(OB(1), OB(2), OB(3), '.k')
% plot3([OB(1) OB(1)+uB(1)], [OB(2) OB(2)+uB(2)], [OB(3) OB(3)+uB(3)],':r')
% plot3([OB(1) OB(1)+vB(1)], [OB(2) OB(2)+vB(2)], [OB(3) OB(3)+vB(3)],':b')
% plot3([OB(1) OB(1)+wB(1)], [OB(2) OB(2)+wB(2)], [OB(3) OB(3)+wB(3)],':g')

end % function defineCoordSystems
