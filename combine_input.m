
% 
% combine_input.m
% 
% White noise input can be used to show that an oscillatory network of LIF
% neurones will tend to respond best when, by chance, there is an
% oscillatory component to the random input. In other words, the
% integration kernel, or transfer function is revealed. This can be seen in
% the spike-triggered average of the input.
% 
% Here, we test the hypothesis that different STA patterns might tell us
% about how two distinct sources of input are combined by a population of
% neurones. Three mechanisms are considered:
% 
%   1) Non-overlapping input. Each input targets a unique set of excitatory
%     neurones.
%   2) Summation. The average of the two inputs is fed to all excitatory
%     neurones.
%   3) Multiplication. The geometric average of the two inputs is fed to
%     all excitatory neurones.
% 
% In addition to varying the mechanism that the network uses to combine
% input, we can also vary what type of input is provided:
%   
%   A) Uniformly distributed white noise. Two independent streams of noise
%     are generated. Each is sampled from a uniform distribution.
%   B) Gaussian distributed white noise. Each input stream is independently
%     sampled from a Gaussian distribution.
%   C) Sinusoidal. Each input is a sine wave. We can systematically vary
%     the relative phase of the two inputs.
%   D) Combination. One input is sinusoidal, the other is white noise.
% 
% Each type of input is tested using each method of combination. Every
% unique pairing of type and method results in a set of STA plots in a 3 by
% 3 arrangement of panels. Each panel plots STAs. Each row of panels show
% STAs for different methods of combination (1-3). Columns show STAs
% computed from different combinations of the input. Col 1 shows STAs
% separately for each input, Col 2 shows the STA of the averaged input, and
% Col 3 shows the STA of the geometric mean of the input. For row 1, non-
% overlapping input, separate STAs are computed for each sub-set of
% excitatory neurones.
% 
% LIF Network parameters aim to replicate those used in Lewis CM, Ni
%   J, Wunderle T, Jendritza P, Lazar A, Diester I, & Fries P. (2020).
%   "Cortical resonance selects coherent input." bioRxiv:
%   2020.2012.2009.417782.
% 
% Written by Jackson Smith - April 2021 - ESI (Fries Lab)
% 


%%% CONSTANTS %%%

% Output directory
OUTDIR = 'C:\Users\smithj\Analysis\Modelling\STA\' ;

% Save figures to PDF file flag, true for save and false for no save.
PDFFLG = true ;

% Duration of input current, in milliseconds
DURONI_MS = 500000 ;

% Half-width of STA, in milliseconds
STAWID_MS = 100 ;

% Input combination method names
INCOMB = { 'Non-overlapping' , 'Summation' , 'Multiplication' } ;


%%% Preparation %%%

% Generate a default network
N = lif( 'network' ) ;

% Point to constants
C = N.C ;

% Convert simulation length from ms to samples
DURONI = ceil( DURONI_MS / C.dt ) ;

% Convert half STA width from ms to samples
STAWID = ceil( STAWID_MS / C.dt ) ;

% Time bin locations in STA at network's sampling rate
timsta = ( -STAWID : +STAWID )' .* C.dt ;

% Allocate matrix of input currents. Each row is a separate time series.
% Columns span time bins.
I = zeros( 2 , DURONI ) ;


%%% Uniform white noise %%%

% Input name
iname = 'Uniform white noise' ;

% Generate two input currents
[ I( 1 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS , 'noise' , [ ] ) ;
[ I( 2 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS , 'noise' , [ ] ) ;

% Evaluate using each method of input combination
sta = runsim( INCOMB , STAWID_MS , iname , N , I ) ;

% Plot and save result
plotsta( OUTDIR , PDFFLG , INCOMB , iname , timsta , sta )


%%% Gaussian white noise %%%

% Input name
iname = 'Gaussian white noise' ;

% Generate two input currents
[ I( 1 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS , 'norm' , [ ] ) ;
[ I( 2 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS , 'norm' , [ ] ) ;

% Evaluate using each method of input combination
sta = runsim( INCOMB , STAWID_MS , iname , N , I ) ;

% Plot and save result
plotsta( OUTDIR , PDFFLG , INCOMB , iname , timsta , sta )


%%% Orthogonal sinusoids %%%

% Input name
iname = 'Sinusoids orthogonal' ;

% Common input parameters, initialise for input 1
par = [ 0 , 1.5 , 40 , 0 ] ;

% Generate input current 1
[ I( 1 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS , 'sine' , par ) ;

% Advance input 2 by pi/2 i.e. 90 degrees
par( 4 ) = pi / 2 ;
[ I( 2 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS , 'sine' , par ) ;

% Evaluate using each method of input combination
sta = runsim( INCOMB , STAWID_MS , iname , N , I ) ;

% Plot and save result
plotsta( OUTDIR , PDFFLG , INCOMB , iname , timsta , sta )


%%% Sinusoid plus uniform white noise %%%

% Input name
iname = 'Sine plus uni noise' ;

% Generate input current 1 - Default sinusoid
[ I( 1 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS ,  'sine' , [ ] ) ;

% Generate input 2 - Uniform white noise
[ I( 2 , : ) , N ] = lif( 'input' , N , 0 , DURONI_MS , 'noise' , [ ] ) ;

% Evaluate using each method of input combination
sta = runsim( INCOMB , STAWID_MS , iname , N , I ) ;

% Plot and save result
plotsta( OUTDIR , PDFFLG , INCOMB , iname , timsta , sta )


%%% DONE %%%

% Clear workspace
clearvars


%%% Script-specific functions %%%

% Compute full set of STAs for each combination method using this network
% and set of input. Returns sta, a 3 x 3 cell array containig STAs. Each
% cell contains a Time x Type array of STA values, with rows spanning time
% bins and columns spanning type of STA. In column 1 of sta, nested arrays
% always have column order [ input 1 STA , input 2 STA ]. All others are
% column vectors averaged across all excitatory neurones, except in row 1,
% where the column order is [ input 1 population , input 2 population ]
% containing STAs computed separately for each set of neurones with
% non-overlapping input (between sets).
function  sta = runsim( INCOMB , w , iname , N , I )
  
  % Report
  fprintf( 'Running simulations for input type: %s\n' , iname )
  
  % Point to network constants
  C = N.C ;
  
  % Allocate sta to accumulate output
  sta = cell( 3 ) ;
  
  % Row counter
  row = 0 ;
  
  % Input combination mechanisms
  for  M = INCOMB , m = M{ 1 } ;  row = row + 1 ;
    
    % Report
    fprintf( '  %s:' , m )
    
    % Set STA flag saying we should group neurones by type
    avgflg = 'type' ;
    
    % Combine input, returns combination in c.
    switch  m
      
      case  'Non-overlapping'
        
        % Compute STA separately for each neurone
        avgflg = 'each' ;
        
        % Allocate a separate input current time series for each neurone in
        % the network
        c = zeros( C.N , size( I , 2 ) ) ;
        
        % Copy input 1 to set 1, and input 2 to set 2
        for  n = 1 : 2 : C.N , c( n , : ) = I( 1 , : ) ; end
        for  n = 2 : 2 : C.N , c( n , : ) = I( 2 , : ) ; end
        
        % I* triggers their computation
        Isum = [ ] ;
        Igeo = [ ] ;
        
      case  'Summation'
        
        % Take the average of the two inputs
        c = mean( I , 1 ) ;
        
        % Point to this for STA, empty Igeo triggers its computation
        Isum =  c ;
        Igeo = [ ] ;
        
      case  'Multiplication'
        
        % Geometric mean of two inputs
        c = sqrt( prod( I , 1 ) ) ;
        
        % Point to this for STA, empty Isum signals its computation
        Isum = [ ] ;
        Igeo =  c ;
      
    end % combo input
    
    % Run simulation
    S = lif( 'sim' , N , c , false ) ;
    
    
    %- Compute STA for each input -%
    
    % Separately for each input
    fprintf( ' Separate' )
    sta{ row , 1 } = [ getsta( N , I( 1 , : ) , S , w , avgflg ) , ...
                       getsta( N , I( 2 , : ) , S , w , avgflg ) ] ;
    
    % Compute summation of inputs?
    if  isempty( Isum ) , Isum = mean( I , 1 ) ; end
    
    % Take STA from average input
    fprintf( ', Average' )
    sta{ row , 2 } = getsta( N , Isum , S , w , avgflg ) ;
      
    % Compute product of inputs?
    if  isempty( Igeo ) , Igeo = sqrt( prod( I , 1 ) ) ; end
    
    % Take STA from geo mean input
    fprintf( ', GeoMean\n' )
    sta{ row , 3 } = getsta( N , Igeo , S , w , avgflg ) ;
    
  end % mechanisms
end % runsim


% Calculate STA across all excitatory units, or separately for each sub-set
function  A = getsta( N , I , S , w , avgflg )
  
  % Compute STAs
  A = lif( 'sta' , N , I , S , w , avgflg ) ;
  
  % Extract STA from A
  switch  avgflg
    
    % Across all excitatory
    case  'type' , A = A( 1 , : )' ;
      
    % Separately for each sub-set of neurones
    case  'each'
      
      A = [ mean( A( N.e( 1 : 2 : end ) , : ) , 1 ) ;
            mean( A( N.e( 2 : 2 : end ) , : ) , 1 ) ]' ;
      
  end % extract
  
end % getsta


% Generates plot of the STAs
function  plotsta( OUTDIR , PDFFLG , INCOMB , iname , t , sta )
  
  % Axes parameters
  AXPARS = { 'XTick' , -60 : 20 : +60 , 'XTickLabel' , [ ] } ;
  
  % Type of STA
  STATYP = { 'Input 1 & 2' , 'Average input' , 'Geo.mean input' } ;
  
  % Find data points to throw away when computing y-axis limits
  nyi = abs( t ) <= 0.5 ;
  
  % Allocate axis limits
  ylims = cell( size( sta ) ) ;
  
  % Rows
  for  r = 1 : size( sta , 1 )
    
    % Concatenate all STAs together
    X = [ sta{ r , : } ] ;
    
    % Discard data points
    X( nyi , : ) = [ ] ;
    
    % Turn into a vector
    X = X( : ) ;
    
    % Determine y-limit based on max and min values
    y = [ min( X ) , max( X ) ] ;
    
    % Expand range by 10%, equally up and down
    y = 1.1 * diff( y ) .* [ -0.5 , +0.5 ]  +  mean( y ) ;
    
    % Copy into accumulator
    ylims( r , : ) = repmat( { y } , 1 , size( sta , 2 ) ) ;
    
  end % y-lims
  
  % Thanks to column-major linear indexing, we need to take the transpose
  % of sta and ylims so that we progress along each row.
    sta = sta' ;
  ylims = ylims' ;
  
  % Determine name of PDF file
  fignam = fullfile( OUTDIR , [ iname , '.pdf' ] ) ;
  
  % Pre-configured figure that is approximately A4
  fig = makfig( -0.95 ) ;
  
  % Panels
  for  pix = 1 : numel( sta )
    
    % New panel
    ax = makax( subplot( 3 , 3 , pix ) , AXPARS{ : } ) ;
    
    % Use default x-axis limit
    ax.XLim = [ -60 , +60 ] ;
    
    % Use measured y-axis limits
    ax.YLim = ylims{ pix } ;
    
    % Plot STAs
    h = plot( ax , t , sta{ pix } , 'k' , 'LineWidth' , 1 ) ;
    
    %- Special cases -%
    
    % Computed for each sub-set of units, and for each set of inputs.
    if  pix == 1
      
      % Make each input a different colour.
      set( h( 3 : 4 ) , 'Color' , ax.ColorOrder( 1 , : ) )
      
      % And then make population 2 a different texture
      set( h( [ 2 , 4 ] ) , 'LineWidth' , 0.5 )
      
    % Computed separately for each subset of neurones
    elseif  pix <= 3
      
      % And then make population 2 a different texture
      set( h( 2 ) , 'LineWidth' , 0.5 )
      
    % Computed for each separate input
    elseif  numel( h ) == 2
      
      % Make each input a different colour.
      set( h( 2 ) , 'Color' , ax.ColorOrder( 1 , : ) )
      
    end % special cases
    
    % Labels
    
    % Top panels tell us what type of STA it is
    if  pix <= 3 , title( ax , STATYP{ pix } ) , end
    
    % Bottom, left panel has axis labels
    if  pix == 2 * 3 + 1
      ax.XTickLabel = ax.XTick ;
      xlabel( 'Time from spike (ms)' )
      ylabel( 'STA (nA)' )
    end
    
    % Not end of row, to next panel
    if  mod( pix , 3 ) , continue , end
    
    % Row number
    r = pix / 3 ;
    
    % Add label that tells us the combination mechanism
    text( ax , ax.XLim( 2 ) , mean( ax.YLim ) , INCOMB{ r } , ...
      'FontSize' , 12 , 'Rotation' , 90 , ...
        'HorizontalAlignment' , 'center' , ...
          'VerticalAlignment' , 'top' )
    
  end % panels
  
  % Switch focus to top-left panel
  ax = subplot( 3 , 3 , 1 ) ;
  
  % Get title object
  h = ax.Title ;
  
  % Bottom location of string
  y = sum( h.Extent( [ 2 , 4 ] ) ) ;
  
  % Print type of input
  text( ax , ax.XLim( 1 ) , y , iname , 'FontSize' , 14 , ...
    'VerticalAlignment' , 'bottom' )
  
  % Show progress
  drawnow
  
  % Print figure to PDF file
  if  PDFFLG
    print( fig , '-bestfit' , '-dpdf' , '-painters' , fignam )
  end
  
end % plotsta

