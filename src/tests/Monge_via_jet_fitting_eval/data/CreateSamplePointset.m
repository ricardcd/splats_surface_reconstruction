% Parameters
OutputFilePath = './x2y2.pts' ;
[X,Y] = meshgrid( -1:.1:1, -1:.1:1 ) ; % Evaluation the grid
Z = -X.^2+-Y.^2 ; % Function to apply
% Z = 3*(1-X).^2.*exp(-(X.^2) - (Y+1).^2) ...
%    - 10*(X/5 - X.^3 - Y.^5).*exp(-X.^2-Y.^2) ...
%    - 1/3*exp(-(X+1).^2 - Y.^2); % peaks
NoiseLevel = 0.5 ;

% RotationAngles = [ 5 5 10 ] ;
% Translation = [ 10 5 -3 ] ;
RotationAngles = [ 0 50 25 ] ;
Translation = [ 0 0 10 ] ;

%% Script

% Display the perfect surface
figure(1) ; clf ;
surf( X, Y, Z, 'FaceAlpha', 0.1 ) ;

% Sort points so first point is the middle one
xs = X(:) ;
ys = Y(:) ;
zs = Z(:) ;

midInd = round( numel( xs ) / 2 ) ;

aux = xs( midInd ) ; xs( midInd ) = xs( 1 ) ; xs( 1 ) = aux ;
aux = ys( midInd ) ; ys( midInd ) = ys( 1 ) ; ys( 1 ) = aux ;
aux = zs( midInd ) ; zs( midInd ) = zs( 1 ) ; zs( 1 ) = aux ;

% Add noise
if NoiseLevel > 0,
  noiseX = NoiseLevel .* rand( numel(xs), 1 ) ;
  noiseY = NoiseLevel .* rand( numel(ys), 1 ) ;
  noiseZ = NoiseLevel .* rand( numel(zs), 1 ) ;
  xs = xs + noiseX ;
  ys = ys + noiseY ;
  zs = zs + noiseZ ;
end


% Apply rotation and translation
xs = xs + Translation(1) ;
ys = ys + Translation(2) ;
zs = zs + Translation(3) ;

RotationMatrix = MobileRPY2RotMat( RotationAngles( 1 ), RotationAngles( 2 ), RotationAngles( 3 ) ) ;
RotatedPts = RotationMatrix * [ xs'; ys'; zs' ] ;

xs = RotatedPts( 1, : )' ;
ys = RotatedPts( 2, : )' ;
zs = RotatedPts( 3, : )' ;

% Save the results on file
fid = fopen( OutputFilePath, 'w' ) ;
for i = 1 : numel( xs ),
   fprintf( fid, '%g %g %g\n', xs(i), ys(i), zs(i) ) ;
end
fclose( fid ) ;

% Show the noisy pointset
% figure(2) ; clf ;
hold on ;
plot3( xs, ys, zs, '.' ) ;
hold on; plot3( xs(1), ys(1), zs(1), 'xr' ) ;

[ minBB maxBB ] = GetAABB( [xs ys zs ] ) ;
