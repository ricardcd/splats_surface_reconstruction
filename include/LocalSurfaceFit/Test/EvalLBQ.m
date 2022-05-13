function z = EvalLBQ( x, y, coeffs ) 

numPts = size( x, 1 ) ;
for i = 1 : numPts,
  z = coeffs( 1 ) * x(i)^2 + coeffs( 2 ) * y(i)^2 + coeffs( 3 ) * x(i) + coeffs( 4 ) * y(i) + coeffs( 5 ) * x(i) * y(i) + coeffs( 6 ) ;  
end




