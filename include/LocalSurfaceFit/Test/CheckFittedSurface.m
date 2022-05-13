%% Parameters
InputLBQFilePath = './LBQ.txt' ;
InputDataPointsFilePath = 'DataPointsCorrupted.xyz' ;
RangeX = [ -1 1 ] ;
RangeY = [ -1 1 ] ;
WeightsH = 2 ; % Change this according to the command line given!
SizeLBQ = 0.5 ;
CircleRadialSteps = 10 ;
CirclePerimeterSteps = 8 ;
OutputLBQFilePath = './DrawnLBQ.obj' ;

%% Load data
LBQ = load( InputLBQFilePath ) ;
if size( LBQ, 2 ) ~= 18,
  error( 'Number of parameters describing an LBQ are 18!' ) ;
end

Data = load( InputDataPointsFilePath ) ;
numDataPoints = size( Data, 1 ) ;

%% Plot results

disp( 'LS Coeficients computed by c++ function: ' ) ;
LBQ( 13:end )' 

% Points
figure ;
hold on ;
plot3( Data( :, 1 ), Data( :, 2 ), Data( :, 3 ), 'xb' ) ;
DrawCircularDiscretizedLBQ( OutputLBQFilePath, LBQ( 1:3 ), LBQ( 4:6 ), LBQ( 7:9 ), LBQ( 10:12 ), LBQ( 13:end ), SizeLBQ, CircleRadialSteps, CirclePerimeterSteps ) 
axis equal ;
title( 'LS result' ) ;
hold off ;

