%% Parameters
InputCoeffsFilePath = './CoeffsRANSAC.txt' ;
InputDataPointsFilePath = 'DataPointsCorrupted.xyz' ;
RangeX = [ -1 1 ] ;
RangeY = [ -1 1 ] ;
WeightsH = 2 ; % Change this according to the command line given!

%% Load data
Coeffs = load( InputCoeffsFilePath ) ;
if size( Coeffs, 2 ) ~= 6,
  error( 'Number of coefficients is 6!' ) ;
end
CoeffsDirect   = Coeffs( 1, : ) ;

Data = load( InputDataPointsFilePath ) ;
numDataPoints = size( Data, 1 ) ;

%% Compute centroid for weighting
if size( Coeffs, 1 ) > 1
  CoeffsWeighted = Coeffs( 2, : ) ;
  EvalPoint = mean( Data( :, 1:2 ) ) ;
  Weights = zeros( numDataPoints, 1 ) ;
  for i = 1 : numDataPoints,
    dist = norm( EvalPoint - Data( i, 1:2 ) ) ;
    Weights( i ) = LSWeightGaussian( dist, WeightsH ) ;
  end
end

%% Evaluate a grid
[ EvalX, EvalY ] = meshgrid( min( RangeX ):.1:max( RangeX ), min( RangeY ):.1:max( RangeY ) ) ;

ZDirect = EvalBivariateQuadratic( CoeffsDirect, [ EvalX(:) EvalY(:) ] ) ;
EvalZDirect = reshape( ZDirect, size( EvalX ) ) ;

if size( Coeffs, 1 ) > 1
  ZWeighted = EvalBivariateQuadratic( CoeffsWeighted, [ EvalX(:) EvalY(:) ] ) ;
  EvalZWeighted = reshape( ZWeighted, size( EvalX ) ) ;
end

%% Plot results

disp( 'LS Coeficients computed by c++ function: ' ) ;
CoeffsDirect' 
disp( 'LS Coeficients computed by matlab function: ' ) ;
CoeffsMat = FitBivariateQuadratic( Data ) 

if size( Coeffs, 1 ) > 1
  disp('')
  disp( ' Centroid = ' ) ;
  EvalPoint 
  disp( 'WLS Coeficients computed by c++ function: ' ) ;
  CoeffsWeighted' 
  disp( 'WLS Coeficients computed by matlab function: ' ) ;
  CoeffsMat = FitWeightedBivariateQuadratic( Data, Weights ) 
end

% Points
figure ;
hold on ;
surf( EvalX, EvalY, EvalZDirect ) ;
plot3( Data( :, 1 ), Data( :, 2 ), Data( :, 3 ), 'xb' ) ;
axis equal ;
title( 'LS result' ) ;
hold off ;

if size( Coeffs, 1 ) > 1
  figure ;
  hold on ;
  surf( EvalX, EvalY, EvalZWeighted ) ;
  plot3( Data( :, 1 ), Data( :, 2 ), Data( :, 3 ), 'xb' ) ;
  axis equal ;
  title( 'WLS result' ) ;
  hold off ;
end

