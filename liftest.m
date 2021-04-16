
% 
% Calculate STA for LIF networks with different amounts of recurrent
% excitatory drive.
% 

% Default parameter set
C = lif( 'default' ) ;

% Replicate to develop different parameter sets
C = repmat( C , 1 , 3 ) ;

% Set 1 is default as stated in Lewis paper, as returned by lif( )

% Set 2 is similar to the example network that Lewis showed me, with no
% recurrent excitatory drive
C( 2 ).V.threshold = -30 ;
C( 2 ).V.reset     = -80 ;
C( 2 ).psp.ee      = +0.0 ;
C( 2 ).psp.ie      = +0.6 ;
C( 2 ).psp.ei      = -0.8 ;
C( 2 ).psp.ii      = -0.8 ;

% Set 3 is the same as set 2, but now we have recurrent excitatory drive
C( 3 ) = C( 2 ) ;
C( 3 ).psp.ee = +0.4 ;

% Millisecond length of stimulus
ms = 5e4 ;

% Width of one frequency bin in discrete Fourier transform
fwid = 1 / ms .* ( 1e3 / C( 1 ).dt ) ;

% Number of bins in 20Hz range, rounded up
fnum = ceil( 20 / fwid ) ;

% Fourier transform smoothing kernel, 5Hz standard deviation
KRNFDT = normpdf( -fnum : +fnum , 0 , 5 )' ;
KRNFDT = KRNFDT ./ sum( KRNFDT ) ;

% CCG convolution kernel from Kohn papers, use for STA
KRNSTA = [ 0.05, 0.25, 0.40, 0.25, 0.05 ]' ;

% Axes parameters
AXPARS = { 'TickDir' , 'out' , 'LineWidth' , 1 , 'FontSize' , 12 , ...
  'Box' , 'off' , 'XLimSpec' , 'tight' , 'XLimMode' , 'auto' , ...
    'YLimSpec' , 'tight' , 'YLimMode' , 'auto' , 'NextPlot' , 'add' } ;

% Screen size
scrsiz = get( groot , 'ScreenSize' ) ;

% Get common input current and rand gen state
N = lif( 'network' ) ;
[ I , N ] = lif( 'input' , N , 0 , ms , 'noise' , [] ) ;

% STA network rng state
rngsta = N.rng ;

% Create figure
fig = figure ;

% Stretch top-to-bottom of screen
fig.OuterPosition( [ 2 , 4 ] ) = scrsiz( [ 2 , 4 ] ) .* [ 1 , 2 / 3 ] ;

% First panel
ax = subplot( numel( C ) , 2 , 1 , AXPARS{ : } ) ;

% Panel index counter
pix = 0 ;

% Parameter sets
for  c = C
  
  % Get network
  N = lif( 'network' , c ) ;
  
  % Set common rng state
  N.rng = rngsta ;
  
  % Run simulation
  S = lif( 'sim' , N , I ) ;
  
  % PSTH [ Excitatory , Inhibitory ]
  X = [ mean( S.spk( N.e , : ) , 1 ) ;
        mean( S.spk( N.i , : ) , 1 ) ]' ;
  
  % Fourier transform of PSTH
  Y = fft( zscore( X ) .* hamming( size( X , 1 ) ) ) ;
  
  % Power spectrum
  P = abs( Y ) .^ 2 ./ size( Y , 1 ) ;
  
  % Frequency vector
  f = ( 0 : size( Y , 1 ) - 1 )' ./ size( Y , 1 ) .* ( 1e3 / c.dt ) ;
  
  % Smooth spectrum
  if  ~ isempty( KRNFDT ) , P = makconv( P , KRNFDT , 's' ) ; end
  
  % New panel
  pix = pix + 1 ;
  ax = subplot( numel( C ) , 2 , pix , AXPARS{ : } ) ;
  
  % Power spectra for excitatory and inhibitory neurones
  plot( f , P )
  
  % Standard axis limits
  xlim( [ 1 , 150 ] )
  ylim( [ 0 , 6 ] )
  
  % Reverse plotting order so that we can see excitatory data better
  ax.Children = flip( ax.Children ) ;
  
  % Show progress
  drawnow
  
  % Allocate spike-triggered average accumulator and spike counter
  STA = zeros( 2 , 200 / c.dt + 1 ) ; num = [ 0 , 0 ] ;
  
  % Neurones
  for  n = 1 : c.N
    
    % Type of neurone
    if  ismember( n , N.e ) , j = 1 ; else , j = 2 ; end
    
    % Locate spikes' linear indices
    spk = find( S.spk( n , : ) ) ;
    
    % Throw away anything near the edges
    spk( spk < 201 | spk > ms / c.dt - 200 ) = [] ;
    
    % Accumulate STA across spikes
    for i = spk
      STA( j , : ) = STA( j , : )  +  I( 1 , i - 200 : i + 200 ) ;
    end
    
    % Count spikes from this neurone
    num( j ) = num( j ) + numel( spk ) ;
  
  end % units
  
  % STA plot
  pix = pix + 1 ;
  ax = subplot( numel( C ) , 2 , pix , AXPARS{ : } ) ;
  
  % Take average across spikes, and smooth slightly
  STA = STA' ./ num ;
  
  % Smooth slightly
  if  ~ isempty( KRNSTA ) , STA = makconv( STA , KRNSTA , 's' ) ; end
  
  % Excitatory and inhibitory STAs
  plot( ( -100 : c.dt : +100 )' , STA , 'LineWidth' , 0.8 )
  
  % Reverse plotting order so that we can see excitatory data better
  ax.Children = flip( ax.Children ) ;
  
  % Standard axis limits
  ylim( [ 0.7 , 0.9 ] )
  
  % Format key parameters into string
  str = sprintf( [ ' V_{threshold} = %d\n V_{reset} = %d\n' , ...
    ' PSP_{EE} = %.1f\n PSP_{IE} = %.1f\n PSP_{EI} = %.1f\n' , ...
      ' PSP_{II} = %.1f' ] , c.V.threshold , c.V.reset , c.psp.ee , ...
        c.psp.ie , c.psp.ei , c.psp.ii ) ;
	
	% Print params
  text( ax.XLim( 1 ) , ax.YLim( 2 ) , str , 'VerticalAlignment' , 'top' )
  
  % Show progress
  drawnow
  
end % psp vals

% Labels on bottom-right
xlabel( 'Time from spike (ms)' )
ylabel( 'STA' )
legend( 'Inhibitory' , 'Excitatory' )

% Labels on bottom-left
subplot( numel( C ) , 2 , numel( C ) * 2 - 1 )
xlabel( 'Frequency (Hz)' )
ylabel( 'PSTH Power' )

% Done!
clearvars

