% Create sample points for the 

%% Parameters
numRandomPoints = 1000 ;
RangeX = [ -1 1 ] ;
RangeY = [ -1 1 ] ;
OutputFilePath = './DataPointsCorrupted.xyz' ;
OutliersPercent = 50 ;
NoiseStd = 0.03 ;

%% Generate random points following the synthetic function
DataX = min( RangeX ) + ( max( RangeX ) - min( RangeX ) ) .* rand( numRandomPoints, 1 ) ;
DataY = min( RangeY ) + ( max( RangeY ) - min( RangeY ) ) .* rand( numRandomPoints, 1 ) ;

%% Synthetic function to generate sample data from
DataZ = -DataX .^ 2 + -DataY .^ 2 ; % Function to apply (change this to generate other points)
Data = [ DataX, DataY, DataZ ] ;

%% Add random noise to the input points
if NoiseStd > 0,
  Noise = NoiseStd * randn( size( Data ) ) ;
  Data = Data + Noise ;
end

Data = AddOutliersToPointsets( Data, [], OutliersPercent, 5 ) ;

%% Save points
Mat2PtsTxt( OutputFilePath, Data ) ;