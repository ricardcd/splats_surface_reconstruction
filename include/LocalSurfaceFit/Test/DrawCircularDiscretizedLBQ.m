function DrawCircularDiscretizedLBQ( OutputFilePath, origins, axisX, axisY, axisZ, coeffs, Sizes, CircleRadialSteps, CirclePerimeterSteps ) 

if nargin < 8,
  CircleRadialSteps = 3 ;
end
if nargin < 9,
  CirclePerimeterSteps = 8 ;
end

% Create discretized unit circle
u = [ 1 0 ] ;
v = [ 0 1 ] ;
RadialStepSize = 1 / CircleRadialSteps ;
PerimeterStepSize = ( 2*pi ) / CirclePerimeterSteps ;

CirclePts = zeros( CircleRadialSteps*CirclePerimeterSteps, 3 ) ;
for i = 1 : CircleRadialSteps,
  radius = RadialStepSize * i ;
  for j = 1 : CirclePerimeterSteps,
    t = PerimeterStepSize * (j-1) ;
    % Evaluation of the parametric circle equation
    pt = radius * cos( t ) * u + radius * sin( t ) * v ;
    CirclePts( ((i-1)*CirclePerimeterSteps)+j, 1:2 ) = pt ;
  end
end  
CirclePts = [ 0 0 0; CirclePts ] ; % Add the 0 point

% --- Debug (Start) ---
% plot( CirclePts( :, 1 ), CirclePts( :, 2 ), 'x' ) ;
% --- Debug  (End)  ---

xs = CirclePts( :, 1 ) ;
ys = CirclePts( :, 2 ) ;
zs = zeros( size( xs ) ) ;
numGridElem = numel( xs ) ;
numLBQ = size( origins, 1 ) ;

% Compute a fixed triangulation
% dt = DelaunayTri( [ xs ys ] ) ;
% FixTri = dt.Triangulation ;
% numFixTri = size( FixTri, 1 ) ;

% First part of the triangulation (triangle fan with the center)
FixTri = zeros( CirclePerimeterSteps + ( (CircleRadialSteps-1)*CirclePerimeterSteps*2 ), 3 ) ;
for j = 2 : CirclePerimeterSteps+1,
  third = j+1 ;
  if third > CirclePerimeterSteps+1,
    third = 2 ;
  end
  FixTri( j-1, : ) = [ 1 j third ] ;
end

curTriInd = CirclePerimeterSteps+1 ;
for i = 2 : CircleRadialSteps,  
  prevCircleIniPt = (CirclePerimeterSteps*(i-2))+2 ; % First point on the last already triangulated circle
  curCircleIniPt = (CirclePerimeterSteps*(i-1))+2 ; % First point on the current circle
  prevCircleCurPt = prevCircleIniPt ;
  curCircleCurPt = curCircleIniPt ;
  for j = 2 : CirclePerimeterSteps+1,
    curCircleNextPt = curCircleCurPt + 1 ;
    prevCircleNextPt = prevCircleCurPt + 1 ;
    if prevCircleNextPt >= ( prevCircleIniPt+CirclePerimeterSteps ),
      prevCircleNextPt = prevCircleIniPt ;
      curCircleNextPt = curCircleIniPt ;
    end
    % Add 2 triangles to the set
    FixTri( curTriInd, : ) = [ prevCircleCurPt curCircleCurPt curCircleNextPt ] ;
    FixTri( curTriInd+1, : ) = [ prevCircleCurPt curCircleNextPt prevCircleNextPt ] ;
    
    % --- Debug (Start) ---
%     triplot( FixTri( 1 : curTriInd+1, : ), CirclePts( :, 1 ), CirclePts( :, 2 ) ) ;
    % --- Debug  (End)  ---
    
    curTriInd = curTriInd + 2 ;
    
    curCircleCurPt = curCircleNextPt ;
    prevCircleCurPt = prevCircleNextPt ;        
  end
end
% --- Debug (Start) ---
% triplot( FixTri, CirclePts( :, 1 ), CirclePts( :, 2 ) ) ;
% --- Debug  (End)  ---

numFixTri = size( FixTri, 1 ) ; % Number of triangles

% Normalize Monges vectors
for i = 1 : numLBQ,
  axisX( i, : ) = axisX( i, : ) ./ norm( axisX( i, : ) ) ;
  axisY( i, : ) = axisY( i, : ) ./ norm( axisY( i, : ) ) ;
  axisZ( i, : ) = axisZ( i, : ) ./ norm( axisZ( i, : ) ) ;
end

% figure ; 
% hold on ;
Pts = zeros( numGridElem*numLBQ, 3 ) ;
Triangles = zeros( numFixTri*numLBQ, 3 ) ;
CurPtsInd = 1 ;
CurTriInd = 1 ;
% Centers = zeros( size( origins ) ) ;
for i = 1 : numLBQ,
    
  % Scale the grid
  xsCur = xs .* ( 2*Sizes( i ) ) ;
  ysCur = ys .* ( 2*Sizes( i ) ) ;
%   xsCur = xs ;
%   ysCur = ys ;
  
  % Create transformation matrix
  Rotation = [ axisX( i, : ) 0 ; ...
               axisY( i, : ) 0 ; ...
               axisZ( i, : ) 0 ; ...
               0 0 0 1 ] ;
  
  Translation = [ eye( 3, 3 ) -origins( i, : )'; 0 0 0 1 ] ;
  
  Transform = Rotation * Translation ;
  
  % Transform points to the current monge coordinates
%   ptsLBQOrigin = Transform * [ xsCur'; ysCur'; zs'; ones( 1, numGridElem ) ] ;
  
  ptsLBQOrigin = [ xsCur'; ysCur'; zs' ] ;

%   figure; 
%   plot3( ptsLBQOrigin( 1, : ), ptsLBQOrigin( 2, : ), ptsLBQOrigin( 3, : ), '.' ) ; 
%   axis equal ; 
%   title( 'ptsLBQOrigin' ) ;
  
  % Evaluate the points on this monge
  ptsLBQ = ptsLBQOrigin( 1:3, : ) ;
  for j = 1 : numGridElem,
    % Substitute Z coordinate by projection on monge
    ptsLBQ( 3, j ) = EvalLBQ( ptsLBQ( 1, j ), ptsLBQ( 2, j ), coeffs ) ;
  end
  
%   figure; 
%   plot3( ptsLBQ( 1, : ), ptsLBQ( 2, : ), ptsLBQ( 3, : ) ) ; 
%   axis equal ; 
%   title( 'ptsLBQ') ;
  
  % Transform back to world coordinates
  ptsWorld = InvHomoPose( Transform ) * [ ptsLBQ ; ones( 1, numGridElem ) ] ;
  
%   origins( i, : )
%   Centers( i, : ) = ptsWorld( 1:3, round( size( ptsWorld, 2 ) / 2 ) )' ;
  
  % Plot points
%   figure;
%   plot3( ptsWorld( 1, : ), ptsWorld( 2, : ), ptsWorld( 3, : ) ) ; 
%   axis equal ; 
%   title( 'ptsWorld' ) ;
  
%   close all ;
  
  Pts( CurPtsInd:CurPtsInd+numGridElem-1, : ) = ptsWorld( 1:3, : )' ;
  Triangles( CurTriInd:CurTriInd+numFixTri-1, : ) = FixTri + ( CurPtsInd -1 ) ;

%   figure ;
%   plot3( Pts( 1:CurPtsInd+numGridElem-1, 1 ), Pts( 1:CurPtsInd+numGridElem-1, 2 ), Pts( 1:CurPtsInd+numGridElem-1, 3 ) ) ; 
%   trimesh( Triangles( 1:CurTriInd+numFixTri-1, : ), Pts( 1:CurTriInd+numFixTri-1, 1 ), Pts( 1:CurTriInd+numFixTri-1, 2 ), Pts( 1:CurTriInd+numFixTri-1, 3 ) ) ; 

  CurPtsInd = CurPtsInd + numGridElem ;
  CurTriInd = CurTriInd + numFixTri ;   
  
end
% Mat2Obj( OutputFilePath, Pts, [], [], Triangles2Faces( Triangles ) ) ;
TriMeshObj( OutputFilePath, Triangles, Pts( :, 1 ), Pts( :, 2 ), Pts( :, 3 ) ) ;
% figure ;
% plot3( Pts( 1, : ), Pts( 2, : ), Pts( 3, : ), '.' ) ; 
% plot3( Centers( :, 1 ), Centers( :, 2 ), Centers( :, 3 ), '.' ) ; 

% --- Debug (Start) ---
%trimesh( Triangles, Pts( :, 1 ), Pts( :, 2 ), Pts( :, 3 ) ) ;
trisurf( Triangles, Pts( :, 1 ), Pts( :, 2 ), Pts( :, 3 ) ) ;
% --- Debug  (End)  ---
