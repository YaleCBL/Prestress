%% ventricle_preStress_main
tic
%%
clear; close all; clc;

%% Plot settings
fontSize=15;
edgeWidth=2;

%% Load
hG2N = 0.000133322;
pressureDiasmmHg = 80;
pressureSysmmHg = 200;
appliedPressureSys= -(pressureSysmmHg*hG2N);

%
%% Material parameter set
%Elastin
c1=0.4; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor
k=c1*k_factor; %Bulk modulus
%Collagen
c1c = 0.5;
c2c = 1;
%
%% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=12; %Optimum number of iterations
max_retries=10; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
%
% Path names
savePath='C:\Users\parik\Desktop\ventricle';
%
% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt'];%Log file name for exporting stress


%% Control parameters

% discQuadMesh parameters
ne=6; %Elements in radius
f=0.5; %Fraction (with respect to outer radius) where central square appears

%Control point parameters
numLayers=18; %Number of control curve layers
heartHeight=65;

%Interpolation settings for warping
interpMethod='linear';
extrapMethod='linear';

%Thickening settings
layerThickness=4; %Wall thickness
numSteps=1; %Number of elements across wall

%% Building a quadrilateral circular template mesh

%Create the mesh
[F_template,V_template]=discQuadMesh(ne,1,f);
V_template(:,3)=0;
F_template=fliplr(F_template);
Eb=patchBoundary(F_template,V_template);
indB=edgeListToCurve(Eb);
indB=indB(1:end-1);

%%
% Smoothen mesh
cPar.n=50;
cPar.Method='LAP';
cPar.RigidConstraints=indB; %Hold on to boundary
[V_template]=patchSmooth(F_template,V_template,[],cPar);

%Define control points on template mesh
V_control_template=V_template(indB,:);

%% Define target control points
% A cell containing a set of control curves is here created. The example
% features a tilted set of ellipse curves that vary in size.

ta=linspace(0,pi/2,numLayers);
a=30+25*(cos(ta)-1); %Radii in first ellipse direction
b=a.*1.2; %Radii in second ellipse direction
p=linspace(0,-heartHeight,numLayers); %Layer offset direction (e.g. Z-coordinate)
% Q=euler2DCM([0.1*pi -0.1*pi 0.25*pi]); %Rotation
Q=euler2DCM(0*[1 1 1]); %Rotation

V_control_target=cell(1,numLayers);
for q=1:1:numLayers
    t = atan2(V_control_template(:,2),V_control_template(:,1));  %Angle
    V_control_target_layer=[a(q)*cos(t) b(q)*sin(t) p(q)*ones(size(t))];
    V_control_target_layer=V_control_target_layer*Q;
    V_control_target{q}=V_control_target_layer;
end

%% Visualize template mesh and control points

cFigure;
subplot(1,2,1); hold on;
title('Template mesh and control points');
hp1=gpatch(F_template,V_template,'w','k',1,edgeWidth);
hp2=plotV(V_control_template,'r.-','LineWidth',3,'MarkerSize',25);
legend([hp1 hp2],{'Template mesh','Template control points'});
axisGeom(gca,fontSize);
camlight headlight;

subplot(1,2,2); hold on;
title('Target control points');

Cp=gjet(numLayers); %Colors for curves
for q=1:1:numLayers
    hp=plotV(V_control_target{q},'k.-','LineWidth',3,'MarkerSize',25);
    hp.Color=Cp(q,:);
end

axisGeom(gca,fontSize);
camlight headlight;
colormap(Cp); caxis([0 numLayers]); icolorbar;
drawnow;

%% Morph template quad mesh

FT=cell(1,numLayers); %Initialize face cell
VT=cell(1,numLayers); %Initialize vertex cell
for q=1:1:numLayers %Loop over layers and process morphing individually
    %Simply copy face description from template
    FT{q}=F_template;

    %Morph vertices
    VT{q}=interpMorph(V_template,V_control_template,V_control_target{q},interpMethod,extrapMethod);
end

%Join face sets (converts cells to normal arrays)
[FT,VT,CT]=joinElementSets(FT,VT);

%%
% Visualizing morphed face sets

cFigure; hold on;
title('Morphed template meshes');
gpatch(FT,VT,CT,'k',1,edgeWidth);

for q=1:1:numLayers
    plotV(V_control_target{q},'g-','LineWidth',3,'MarkerSize',25);
end

axisGeom(gca,fontSize);
camlight headlight;
colormap gjet;  icolorbar;
drawnow;

%% Build hexahedral elements of interior
% Loop over face layers and use faces as tops/bottoms of volumetric
% elements.

E1=[];
C1=[];
for q=1:1:numLayers-1
    e=[FT(CT==q,:) FT(CT==q+1,:)];
    E1=[E1; e];
    C1=[C1; q*ones(size(e,1),1)];
end

[F1,CF1]=element2patch(E1,C1);
V1=VT;

%%
% Visualize interior mesh

cFigure; hold on;
title('Hex mesh fluid');

gpatch(F1,V1,CF1,'k',1,edgeWidth);
% patchNormPlot(F1,V1);

axisGeom(gca,fontSize);
camlight headlight;
colormap gjet;  icolorbar;
drawnow;

%% Bound wall mesh from by offsetting interior

%%
% Get boundary faces

indBoundaryFaces=tesBoundary(F1,V1);
Fb1=F1(indBoundaryFaces,:);
Cb1=CF1(indBoundaryFaces,:);

Cb1_V=faceToVertexMeasure(Fb1,V1,Cb1);

logicSides=~all(Cb1_V(Fb1)==1,2);

Fb1_sides=Fb1(logicSides,:);
Cb1_sides=Cb1(logicSides,:);

%%
% Smooth outer mesh
indSmooth=unique(Fb1_sides(:));
[Ft,Vt]=patchCleanUnused(Fb1_sides,V1);

cPar.n=25;
cPar.Method='HC';
cPar.RigidConstraints=indB;
[Vt]=patchSmooth(Ft,Vt,[],cPar);
V1(Fb1_sides(:),:)=Vt(Ft(:),:);

%%
% Thicken mesh to form layer

[Ft,Vt]=patchCleanUnused(Fb1_sides,V1);
[E2,V2,F2_1,F2_2]=patchThick(Ft,Vt,1,layerThickness,numSteps);
C2=repmat(Cb1_sides,[numSteps,1]);
[F2,CF2]=element2patch(E2,C2);

%%
% Visualize wall mesh
indInner = unique(F2_1);
indBoundFix = find(V2(:,3)>-2.5 );
%
cFigure; hold on;
title('Hex mesh wall');
gpatch(F2,V2,CF2,'k',1,edgeWidth);
plotV(V2(indInner,:),'r.','markersize',20)
plotV(V2(indBoundFix,:),'g.','markersize',20)
patchNormPlot(F2,V2);
axisGeom(gca,fontSize);
camlight headlight;
colormap gjet;  icolorbar;
drawnow;

%% Renaming elements
VT = V2;
F_pressure = F2_1;
ET = E2;
bcSupportList = indBoundFix;
bcPrescribeList = indInner;

%% Defining the FEBio input structure
appliedPressure = appliedPressureSys;

%
for q = 1:1:size(appliedPressure,2)

    %Get a template with default settings
    [febio_spec]=febioStructTemplate;

    %febio_spec version
    febio_spec.ATTR.version='3.0';

    %Module section
    febio_spec.Module.ATTR.type='solid';

    %Control section
    febio_spec.Control.analysis='STATIC';
    febio_spec.Control.time_steps=numTimeSteps;
    febio_spec.Control.step_size=1/numTimeSteps;
    febio_spec.Control.solver.max_refs=max_refs;
    febio_spec.Control.solver.max_ups=max_ups;
    febio_spec.Control.time_stepper.dtmin=dtmin;
    febio_spec.Control.time_stepper.dtmax=dtmax;
    febio_spec.Control.time_stepper.max_retries=max_retries;
    febio_spec.Control.time_stepper.opt_iter=opt_iter;

    %Material section
    materialName1='Material1';
    febio_spec.Material.material{1}.ATTR.name=materialName1;
    febio_spec.Material.material{1}.ATTR.type='Ogden';
    febio_spec.Material.material{1}.ATTR.id=1;
    febio_spec.Material.material{1}.c1=c1;
    febio_spec.Material.material{1}.m1=m1;
    febio_spec.Material.material{1}.c2=c1;
    febio_spec.Material.material{1}.m2=-m1;
    febio_spec.Material.material{1}.k=k;

    % Mesh section
    % -> Nodes
    febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(VT,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=VT; %The nodel coordinates

    % -> Elements
    partName1='Part1';
    febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
    febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
    febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(ET,1))'; %Element id's
    febio_spec.Mesh.Elements{1}.elem.VAL=ET; %The element matrix

    % -> Surfaces
    surfaceName1='LoadedSurface';
    febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
    febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
    febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;

    % -> NodeSets

    nodeSetName1='bcSupportList';
    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

    %MeshDomains section
    febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
    febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

    %Boundary condition section
    % -> Fix boundary conditions
    febio_spec.Boundary.bc{1}.ATTR.type='fix';
    febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{1}.dofs='x,y,z';


    %Loads section
    % -> Surface load
    febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
    febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
    febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
    febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;
    febio_spec.Loads.surface_load{1}.pressure.VAL=1;
    %
    %LoadData section
    % -> load_controller
    febio_spec.LoadData.load_controller{1}.ATTR.id=1;
    febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
    febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
    febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 appliedPressure(q)];

    %Output section
    % -> log file
    febio_spec.Output.logfile.ATTR.file=febioLogFileName;
    febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
    febio_spec.Output.logfile.node_data{1}.VAL=1:size(VT,1);

    febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
    febio_spec.Output.logfile.element_data{1}.ATTR.data='sx;sy;sz;sxy;syz;sxz';
    febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
    febio_spec.Output.logfile.element_data{1}.VAL=1:size(ET,1);
    %
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    %
    %% Running the FEBio analysis

    febioAnalysis.run_filename=febioFebFileName; %The input file name
    febioAnalysis.run_logname=febioLogFileName; %The name for the log file
    febioAnalysis.disp_on=1; %Display information on the command window
    febioAnalysis.runMode='external';
    febioAnalysis.maxLogCheckTime=10; %Max log file checking time

    [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

    if runFlag==1 %i.e. a succesful run
        %Importing displacement
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
        N_disp_mat=dataStruct.data; %Displacement
        U{q} = N_disp_mat(:,:,end);
        %
        %import elemental stresses from logfile
        dataStruct_stress = importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
        E_stress_disp = dataStruct_stress.data; %Elemental stresses
        E_stress_mat{q}= E_stress_disp(:,:,end);

    end
    %Creating systolic geometries
    %
    V_DEF_int =  VT + U{1,q};
    V_DEF_sys{q} = V_DEF_int;
    %Inner layer of systolic geometry for virtual work
    V_sysInner = V_DEF_int(bcPrescribeList,:);
    %
end


%% Defining the deformed configuration as the configuration from medical images and providing a perturbation
V_new_ref = VT+U{1,end};%
%
norm_U = sqrt(sum((U{1,end}).^2, 2));
%
%Creating a perturbed geometry
U_perturb = ((U{1,end})./norm_U);
rowsWithNaN = any(isnan(U_perturb), 2);
% % Replace those rows with ones
U_perturb(rowsWithNaN, :) = 0;

%%
V_perturbed =  V_new_ref - 0.0001*U_perturb;
%Inner layer of systolic geometry for virtual work
V_sysInner_perturbed = V_perturbed(bcPrescribeList,:);
%

%%
cFigure; hold on;
gpatch(F2,V_new_ref,'none','k',0.5);
gpatch(F2,V_perturbed,'kw','r',0.5);
axisGeom
%%
VT_img = V_new_ref;
%%
%% FIBER ORIENTATION ASSIGNMENT: Circumferential fibers
%% Step 1: Define Gauss points for 2x2x2 integration
a = 1/sqrt(3);
gp = [-a -a -a;
       a -a -a;
       a  a -a;
      -a  a -a;
      -a -a  a;
       a -a  a;
       a  a  a;
      -a  a  a]; % (ξ,η,ζ) coordinates of 8 Gauss points

numElements = size(ET,1);
numGP = size(gp,1);

fiber_dir_IP = zeros(numElements,numGP,3);
pos_IP = zeros(numElements,numGP,3);

%% Step 2: Loop over elements and compute integration point coordinates
for e = 1:numElements
    % Node coordinates of the element
    Xe = VT_img(ET(e,:),:); % 8x3
    
    for g = 1:numGP
        xi = gp(g,1); eta = gp(g,2); zeta = gp(g,3);
        
        % Shape functions for linear hex8
        N = 1/8 * [(1-xi)*(1-eta)*(1-zeta);
                   (1+xi)*(1-eta)*(1-zeta);
                   (1+xi)*(1+eta)*(1-zeta);
                   (1-xi)*(1+eta)*(1-zeta);
                   (1-xi)*(1-eta)*(1+zeta);
                   (1+xi)*(1-eta)*(1+zeta);
                   (1+xi)*(1+eta)*(1+zeta);
                   (1-xi)*(1+eta)*(1+zeta)];
        
        % Global position of integration point
        x_IP = N' * Xe;
        pos_IP(e,g,:) = x_IP;
        
        % Compute circumferential direction at that point
        xc = x_IP(1); yc = x_IP(2);
        v_radial = [xc, yc, 0];
        
        if norm(v_radial) < 1e-8
            v_circ = [1,0,0]; % default
        else
            z_axis = [0,0,1];
            v_circ = cross(z_axis, v_radial);
            v_circ = v_circ ./ norm(v_circ);
        end
        
        fiber_dir_IP(e,g,:) = v_circ;
    end
end

%% Step 3: (Optional) Visualize fibers at integration points
cFigure; hold on;
title('Circumferential fibers at integration points');
gpatch(F2,VT_img,'w','none',0.1);
quiver3(pos_IP(:,:,1),pos_IP(:,:,2),pos_IP(:,:,3),...
        fiber_dir_IP(:,:,1),fiber_dir_IP(:,:,2),fiber_dir_IP(:,:,3),...
        0.3,'r','LineWidth',1);
axisGeom;
camlight headlight;
drawnow;

%% Step 4: Save to struct
fiberField.hexIP_dir = fiber_dir_IP;
fiberField.hexIP_pos = pos_IP;
save('fiberField_hex8_IP_circumferential.mat','fiberField');

%% COMPUTE CAUCHY STRESS AT INTEGRATION POINTS (Fiber-reinforced model)

numElements = size(fiberField.hexIP_dir,1);
numGP = size(fiberField.hexIP_dir,2);

% Material parameters
lc1 = 1.062;
c1c = 0.5;
c2c = 1.0;

% Initialize storage
stress_Cauchy = zeros(numElements,numGP,3,3);  % Cauchy stress tensor at each IP
I1_c1_all = zeros(numElements,numGP);
I4_c1_all = zeros(numElements,numGP);

%% Loop over all integration points
for e = 1:numElements
    for g = 1:numGP
        % Fiber direction at integration point
        mc1 = squeeze(fiberField.hexIP_dir(e,g,:))'; % (1x3)
        mc1 = mc1 / norm(mc1); % Ensure unit vector
        
        % Compute Gc_c1 (fiber stretch tensor)
        Gc_c1 = lc1*(mc1'*mc1) + (1/sqrt(lc1))*(eye(3) - (mc1'*mc1));
        
        % Right-Cauchy deformation tensor
        Cc_c1 = Gc_c1' * Gc_c1;
        
        % Invariants
        I1_c1 = trace(Cc_c1);
        I4_c1 = dot(mc1',(Cc_c1*mc1'));
        
        % Store invariants
        I1_c1_all(e,g) = I1_c1;
        I4_c1_all(e,g) = I4_c1;
        
        % Second Piola-Kirchhoff (fiber) stress
        S_col_c1 = c1c * (I4_c1 - 1) * exp( c2c * (I4_c1 - 1)^2 ) * (mc1'*mc1);
        
        % Cauchy stress tensor
        sigma_c1 = Gc_c1 * S_col_c1 * Gc_c1';
        
        stress_Cauchy(e,g,:,:) = sigma_c1;
    end
end

%Store results
stressResults.Cauchy = stress_Cauchy;
%%
%% INTERNAL FORCE COMPUTATION AND VISUALIZATION
% Computes F_int = ∫ σ ∇N dV at each node (for HEX8 elements)
% Requires:
%   VT_img          : nodal coordinates [nNodes x 3]
%   ET              : element connectivity [nElem x 8]
%   stress_Cauchy   : [nElem x 8 x 3 x 3] Cauchy stress at Gauss points

clearvars -except VT_img ET stress_Cauchy F2 F_pressure appliedPressureSys savePath
clc

numElements = size(ET,1);
numNodes = size(VT_img,1);

%% --- 2x2x2 Gauss integration points and weights ---
a = 1/sqrt(3);
gp = [-a -a -a;
       a -a -a;
       a  a -a;
      -a  a -a;
      -a -a  a;
       a -a  a;
       a  a  a;
      -a  a  a];
w = ones(8,1);

%% --- Initialize global internal force vector ---
Fint_global = zeros(numNodes,3);

%% --- Loop over all elements ---
for e = 1:numElements
    Xe = VT_img(ET(e,:),:); % element nodal coordinates (8x3)
    Fint_e = zeros(8,3);     % element internal force vector
    
    for g = 1:8
        xi = gp(g,1); eta = gp(g,2); zeta = gp(g,3);
        
        % Derivatives of shape functions wrt local coords (8x3)
        dN_dxi = 1/8 * [
            -(1-eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1-eta);
             (1-eta)*(1-zeta),  -(1+xi)*(1-zeta),  -(1+xi)*(1-eta);
             (1+eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1+eta);
            -(1+eta)*(1-zeta),   (1-xi)*(1-zeta),  -(1-xi)*(1+eta);
            -(1-eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1-eta);
             (1-eta)*(1+zeta),  -(1+xi)*(1+zeta),   (1+xi)*(1-eta);
             (1+eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1+eta);
            -(1+eta)*(1+zeta),   (1-xi)*(1+zeta),   (1-xi)*(1+eta)];
        
        % Jacobian matrix and determinant
        J = Xe' * dN_dxi;
        detJ = det(J);
        if detJ < 0
            warning('Negative Jacobian in element %d',e);
        end
        
        % Gradient of shape functions wrt global coords
        dN_dX = dN_dxi / J;  % (8x3)
        
        % Stress tensor at current integration point
        sigma_gp = squeeze(stress_Cauchy(e,g,:,:)); % (3x3)
        
        % Integration weight
        wgp = w(g) * detJ;
        
        % Compute nodal force contributions
        for a = 1:8
            gradNa = dN_dX(a,:)'; % (3x1)
            Fint_e(a,:) = Fint_e(a,:) + (sigma_gp * gradNa)' * wgp;
        end
    end
    
    % --- Assemble element forces into global vector ---
    for a = 1:8
        nodeID = ET(e,a);
        Fint_global(nodeID,:) = Fint_global(nodeID,:) + Fint_e(a,:);
    end
end

%% --- Visualization of Internal Force Vectors ---
Fmag = sqrt(sum(Fint_global.^2,2));
scaleFactor = 1;  % adjust this for clarity in quiver3 arrows

cFigure; hold on;
title('Internal Nodal Force Vectors (∫σ∇N dV)');
if exist('F2','var')
    gpatch(F2,VT_img,'w','none',0.3); % ventricular surface if available
else
    plotV(VT_img,'k.'); % fallback if surface not defined
end
quiver3(VT_img(:,1), VT_img(:,2), VT_img(:,3), ...
        -Fint_global(:,1), -Fint_global(:,2), -Fint_global(:,3), ...
        scaleFactor, 'r', 'LineWidth',1.5);

colormap turbo;
colorbar;
axisGeom;
camlight headlight;
drawnow;

%% --- Optional: Save results ---
internalWorkResults.Fint_global = Fint_global;
internalWorkResults.Fint_magnitude = Fmag;
save('internalForceVectors.mat','internalWorkResults');

%% ------------------------------------------------------------------------
% Finite Element Pressure Force Calculation on Inner Ventricular Surface
% -------------------------------------------------------------------------
% Computes equivalent nodal forces due to applied pressure using FE surface
% integration (Gauss quadrature on each 4-node quad).
% Inputs: VT_img (nodal coords), F_pressure (inner faces), appliedPressureSys
% -------------------------------------------------------------------------

fprintf('\n--- Computing finite element pressure forces (quadrature integration) ---\n');

% Safety checks
if ~exist('VT_img','var') || isempty(VT_img)
    error('VT_img not found.');
end
if ~exist('F_pressure','var') || isempty(F_pressure)
    error('F_pressure (inner surface faces) missing.');
end
if ~exist('appliedPressureSys','var') || isempty(appliedPressureSys)
    error('appliedPressureSys (pressure value) missing.');
end

numNodes = size(VT_img,1);
numFaces = size(F_pressure,1);

% Initialize global vector
Fext_pressure = zeros(numNodes,3);

% Pressure value
p_scalar = appliedPressureSys;

% 2x2 Gauss quadrature points and weights
gp = [-1/sqrt(3),  1/sqrt(3)];
w  = [1, 1];

% Loop over all quad faces
for f = 1:numFaces
    nodeIDs = F_pressure(f,:);
    Xe = VT_img(nodeIDs,:); % 4x3 coordinates of face

    % Element nodal force accumulator
    fe = zeros(4,3);

    % Integrate over 2x2 Gauss points
    for i = 1:2
        for j = 1:2
            xi  = gp(i);
            eta = gp(j);
            wt  = w(i)*w(j);

            % Shape functions (bilinear quad)
            N = 0.25 * [(1 - xi)*(1 - eta);
                        (1 + xi)*(1 - eta);
                        (1 + xi)*(1 + eta);
                        (1 - xi)*(1 + eta)];  % 4x1

            % Derivatives w.r.t parent coordinates
            dN_dxi = 0.25 * [-(1 - eta),  (1 - eta),  (1 + eta), -(1 + eta)];
            dN_deta = 0.25 * [-(1 - xi), -(1 + xi),  (1 + xi),  (1 - xi)];

            % Tangent vectors in global coordinates
            dx_dxi  = dN_dxi  * Xe;  % 1x3
            dx_deta = dN_deta * Xe;  % 1x3

            % Surface normal (not normalized)
            normalVec = cross(dx_dxi, dx_deta);
            dS = norm(normalVec);
            if dS == 0, continue; end

            n = normalVec / dS;   % unit normal
            J_s = dS;             % surface Jacobian magnitude

            % Pressure traction vector
            traction = p_scalar * n;  % 1x3

            % Nodal force contribution
            for a = 1:4
                fe(a,:) = fe(a,:) + N(a) * traction * J_s * wt;
            end
        end
    end

    % Assemble element nodal forces to global
    for a = 1:4
        Fext_pressure(nodeIDs(a),:) = Fext_pressure(nodeIDs(a),:) + fe(a,:);
    end
end


%% Visualization
Fpress_mag = sqrt(sum(Fext_pressure.^2,2));
% scaleFactor = 0.2 * max(range(VT_img)) / max(Fpress_mag);

% Compute geometric range manually (no toolbox)
geomRange = max(max(VT_img,[],1) - min(VT_img,[],1)); % max dimension span
scaleFactor = 0.2 * geomRange / max(Fpress_mag);


cFigure; hold on;
title('Finite Element Assembled Nodal Pressure Forces');
if exist('F2','var')
    gpatch(F2,VT_img,'w','k',0.2);
else
    plotV(VT_img,'k.','MarkerSize',4);
end
quiver3(VT_img(:,1),VT_img(:,2),VT_img(:,3), ...
        Fext_pressure(:,1),Fext_pressure(:,2),Fext_pressure(:,3), ...
        0.07*scaleFactor, 'b', 'LineWidth', 1.2);
scatter3(VT_img(:,1),VT_img(:,2),VT_img(:,3),20,Fpress_mag,'filled');
axisGeom; camlight headlight; drawnow;
%%
%% ------------------------------------------------------------------------
% Compute and visualize net nodal forces
% -------------------------------------------------------------------------

fprintf('\n--- Computing net nodal forces ---\n');

% Ensure both force fields exist and have the same size
if ~exist('Fint_global','var') || ~exist('Fext_pressure','var')
    error('Missing Fint_global or Fext_pressure in workspace.');
end

% --- Compute net force per node ---
Fnet_global = -Fext_pressure - Fint_global;  % subtract internal (resisting) from external (driving)
Fnet_mag = sqrt(sum(Fnet_global.^2,2));

% --- Visualization ---
cFigure; hold on;
title('Net Nodal Force Vectors (Pressure + Internal)');
if exist('F2','var')
    gpatch(F2,VT_img,'w','k',0.15);
else
    plotV(VT_img,'k.','MarkerSize',4);
end

% Scale arrows for visualization
geomRange = max(max(VT_img,[],1) - min(VT_img,[],1)); % geometric span
scaleFactor = 0.2 * geomRange / max(Fnet_mag);

% Plot quiver arrows showing direction and magnitude
quiver3(VT_img(:,1), VT_img(:,2), VT_img(:,3), ...
        Fnet_global(:,1), Fnet_global(:,2), Fnet_global(:,3), ...
        0.07*scaleFactor, 'm', 'LineWidth', 1.3);

% Add color-coded magnitude
scatter3(VT_img(:,1),VT_img(:,2),VT_img(:,3),20,Fnet_mag,'filled');
colormap turbo; colorbar;
axisGeom; camlight headlight;
drawnow;
%%
% Inside main_script.m
disp('Running febio_run.m');
run('febio_run_vent.m');  % Runs another_script.m
disp('END OF SCRIPT');
%%




























































