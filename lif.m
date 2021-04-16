
function  varargout = lif( fstr , varargin )
% 
% lif( <function name> , args ... )
% 
% Create and use a network of leaky integrate and fire neurones. The
% function name is one of the following strings:
% 
% 
% C = lif( 'default' )
%   
%   Returns struct C with a default set of network parameters. These are
%   chosen to match, as mutch as possible, values taken from Lewis CM, Ni
%   J, Wunderle T, Jendritza P, Lazar A, Diester I, & Fries P. (2020).
%   "Cortical resonance selects coherent input." bioRxiv:
%   2020.2012.2009.417782. The 'C' is for 'C'onstants, as opposed to
%   variables.
% 
% 
% N = lif( 'network' )
% N = lif( 'network' , C )
%   
%   Returns struct N with the instantiation of a LIF neural network using
%   parameters in struct C. If C is not provided then default parameters
%   are used.
% 
% 
%         I = lif( 'input' , N , off , on )
% [ I , N ] = lif( 'input' , N , off , on , type , par )
% [ I , N ] = lif( 'input' , N , off , on , ... , 'repeat' , rep )
% 
%   Returns array I containing input current in nA for the LIF neurones in.
%   N. I is either a row vector, or a matrix with Neurones indexed across
%   rows. In each case, a single row contains a time series of the input
%   current in chronological order. Hence columns of I index time steps,
%   with C.dt duration each, where C = N.C.
%   
%   Input arguments 'off' and 'on' give durations in milliseconds. The
%   input current contains three epochs [ pre-off , ON , post-off ]. The
%   pre- and post-off epochs are each 'off' milliseconds long; the single
%   ON epoch is 'on' milliseconds long. Durations are converted to time
%   steps through division by C.dt and then rounding up e.g. off_steps =
%   ceil(off/C.dt). Thus, I will have ( 2 * off_steps + on_steps ) / C.dt
%   columns. 'off' can be 0ms, but 'on' must have a positive value.
%   
%   The type of signal that is presented during the ON epoch can be set
%   with the name/value pair given in type and par, where type is a string
%   naming what kind of signal to use, and par is a numeric vector
%   containing parameters. If par is empty i.e. [] then default parameters
%   are used based on the parameters of N. If type and par are ommitted
%   then a default signal is generated.
%   
%   type strings:
%   
%   'const' (default) - A constant value is presented during the ON phase.
%     par is a scalar number giving the constant values, in nA.
%     If par is empty, then 1.5nA is used by default.
%   
%   'ramp' - Current increases (or decreases) linearly from the start to
%     the end of the ON phase. par is [ base , amp ], in nA. If par is
%     empty then base = 0nA and amp = 1.5nA. If amp > 0, then the waveform
%     at time t milliseconds from the start of the ON phase will be:
%       base  +  t / on * amp.
%     If amp < 0 then the waveform will be:
%       base  +  abs( amp )  +  t / on * amp.
%     Thus, if amp > 0 then the ramp is increasing over time, and if amp <
%     0 then the ramp is decreasing over time.
%   
%   'sine' - Sine wave, where par is [ base , amp , freq ]. The baseline
%     (base) and amplitude (amp) are in nA, and the frequency (freq) is in
%     Hz. If par is empty then base = 0nA, amp = 1.5nA, and freq = 40Hz by
%     default. The waveform is constructed so that the current at time t in
%     milliseconds from the start of the ON phase is:
%       base  +  ( 1 + sin( 2*pi * freq * t/1e3 - pi/2 ) ) / 2 * amp.
%     The -pi/2 term causes the sine wave to start at its minimum, avoiding
%     any abrupt step at the start of the ON phase, if base = 0.
%   
%   'noise' - White noise stimulus sampled from a uniform distribution. par
%     is [ base , amp ] in nA so that the current at time point t is
%     sampled from the unif( base , amp ) distribution; in other words, all
%     values lie between base and amp. If par is empty then base = 0nA and
%     amp = 1.5nA.
%     
%   'norm' - White noise stimulus sampled from a normal distrubition i.e. a
%     Gaussian distribution. par is [ avg , sd ] in nA where avg is the
%     average or mean of the distribution, and sd is the standard
%     deviation. If par is empty then avg = 1.5/2 and sd = avg/3.5. Thus,
%     the current at time t is sampled from distribution N( avg , sd ^ 2 ).
%     The sampled values will not be perfectly normal in their distribution
%     because values will be clipped at plus and minus 3.5*sd. In other
%     words, if I( t ) is the current at time t and I( t ) > avg+3.5*sd
%     when first sampled, then the value is replaced with I( t ) =
%     avg+3.5*sd; alternatively, if I( t ) < avg-3.5*sd when sampled then
%     I( t ) = avg-3.5*sd is the replacement.
%   
%   Optional name/value pair 'repeat' and rep can be supplied. rep is a
%   scalar,signed numeric value. If rep = 0 (default) then I will be a row
%   vector, returning a single copy of the input current time series. On
%   the other hand, if rep is non-zero and positive e.g. rep = +1 then the
%   input current will be copied once for each LIF neurone in N, returning
%   I as a matrix with C.N identical rows. A special case is when rep is
%   non-zero and negative e.g. rep = -1. This triggers a special behaviour
%   when using type = 'noise' or 'norm'. In this case, a unique white-noise
%   time series is sampled in the ON phase for each LIF neurone; no two
%   neurones will receive the same input current time series, in this case.
%   
%   Note that white-noise generation calls rand or randn. The random-number
%   generator state in input argument N will be applied before sampling
%   numbers. After sampling numbers, the updated rng state is returned in
%   output argument N.
%   
% 
% [ S , N ] = lif( 'sim' , N , I )
%       ... = lif( 'sim' , N , I , volflg )
%   
%   Runs a simulation using LIF network struct N and input current I in nA
%   for all neurones. I can be a single row vector containing the input
%   current time series, or a matrix with a row of current input values for
%   each separate LIF neurone in N. If I is a row vector then the same
%   input is applied to each neurone. Either way, the input is first scaled
%   by terms C.Ie and C.Ii, according to the type of each neurone.
%   
%   Returns struct S containing the results of the simulation, with fields:
%   
%     S.spk - Neurones x Time logical array. Each row is a spike raster,
%       where 1's mark out the time bins in which a spike was fired by that
%       neurone.
%     
%     S.vol - Neurones x Time double array of membrane voltages.
% 
%   Note that the initial membrane potential will be randomly sampled
%   uniformly between C.V.leak and C.V.threshold. Hence, state N.rng is
%   updated and returned in output N.
%   
%   Optional input argument volflg is used to control whether or not
%   membrane voltages are saved and returned in S.vol. If volflg is true or
%   non-zero (default) then voltages are saved and returned. Otherwise, if
%   volflg is false or zero, then S.vol returns an empty array i.e. [ ].
%   This behaviour may be desirable if voltage membrane output is not
%   required or if the simulation is very large; substantially less memory
%   will be used if the voltages are not saved and returned.
% 
% 
% A = lif( 'sta' , N , I , S )
% A = lif( 'sta' , N , I , S , w )
% A = lif( 'sta' , N , I , S , w , avgflg )
%   
%   Computes spike-triggered average A of input I aligned to spikes in
%   raster S.spk from network N. If I is a vector, then all spikes are
%   compared against the same input. But if I is a matrix then I must have
%   the same size as S.spk, and spikes from row S.spk( r , : ) are only
%   compared against input from I( r , : ). w is optional and gives the
%   width of the STA in milliseconds (default 200ms). The STA is computed
%   from -w up to +w in C.dt steps. Hence, A is ceil( 2 * w / C.dt ) + 1
%   long in the time domain; the +1 term accounts for the time bin that
%   contains the spike. The size of A depends on optional input avgflg (set
%   w to empty i.e. [ ] for default value), which is a character string
%   that controls how the STA is averaged across neurones in N. avgflg can
%   be one of the following:
%   
%      'all' - STA computed from all spikes. Returns row vector in A.
%     'type' - STA computed separately for excitatory and inibitory
%       neurones. A is a 2 x Time array, with row order [ excitatory STA ;
%       inhibitory STA ]. [DEFAULT]
%     'each' - STA computed separately for each neurone. A is C.N x Time
%       array in which A( i , : ) is the STA for the nth neurone in N,
%       corresponding to spike raster S.spk( i , : ).
% 
% 
% Written by Jackson Smith - ESI Fries Lab - April 2021
% 


%%% Check Input %%%

% fstr must be a char vector
if  ~ ischar( fstr )  ||  ~ isvector( fstr )  ||  ~ isrow( fstr )
  error( 'fstr must be a string i.e. char row vector' )
end


%%% Select function %%%

switch  fstr
  
  
  %-- Define a default set of network parameters --%
  
  case  'default'
    
    % Check max number of input/output args
     narginchk( 1 , 1 )
    nargoutchk( 0 , 1 )
    
    % Number of excitatory neurones
    C.Ne = 200 ;
    
    % Number of inhibitory neurones
    C.Ni =  50 ;
    
    % Total number of neurones
    C.N = C.Ne + C.Ni ;
    
    % Scaling term that is applied to input current, separately for
    % excitatory and inhibitory neurones. That is, raw input current I is
    % scaled by C.Ie .* I and C.Ii * I before being delivered to the LIF
    % neurones of either type.
    C.Ie = 1.0 ;
    C.Ii = 0.0 ;
    
    % Maximum post-synaptic membrane potential in a downstream neurone that
    % is induced by an incoming spike from an afferent neurone. The coding
    % <downstream><afferent> gives the type of each neurone, where type is
    % encoded as e - excitatory or i - inhibitory. Units in mV.
    C.psp.ee =  0 ;
    C.psp.ie = +1 ;
    C.psp.ei = -1 ;
    C.psp.ii = -1 ;
    
    % Duration of a single time-step, in milliseconds
    C.dt = 0.5 ;
    
    % Cellular membrane capacitance, in nF
    C.C = 0.5 ;
    
    % Cellular membrane resistance, in Mohms
    C.R = 40 ;
    
    % Membrane time constant
    C.tau = C.R * C.C ;
    
    % Membrane voltage during a spike, in mV
    C.V.spike = +30 ;
    
    % Voltage threshold for triggering a spike, in mV
    C.V.threshold = -40 ;
    
    % Reset voltage immediately after firing a spike, in mV
    C.V.reset = -70 ;
    
    % The equilibrium membrane potential in absence of input, in mV
    C.V.leak = -60 ;
    
    % Store current state of random number generator
    rngtmp = rng ;
    
    % Re-set default random number generator state
    rng( 'default' )
    
    % Store default RNG state
    C.rng = rng ;
    
    % Restore state of rng
    rng( rngtmp )
    
    % Return param struct
    varargout = { C } ;
    
    
  %-- Instantiate a LIF network --%
  
  case  'network'
    
    % Check max number of input/output args
     narginchk( 1 , 2 )
    nargoutchk( 0 , 1 )
    
    % Parameters given explicitly, use them
    if  nargin == 2
      
      C = varargin{ 1 } ;
      
    % No parameters given, fetch default
    else
      
      C = lif( 'default' ) ;
      
    end % params
    
    % Store current state of random number generator
    rngtmp = rng ;
    
    % Set rng state from param set
    rng( C.rng )
    
    % Store network parameters
    N.C = C ;
    
    % Index vector of excitatory neurones
    N.e = 1 : C.Ne ;
    
    % Index vector of inhibitory neurones
    N.i = C.Ne + 1 : C.N ;
    
    % Allocate weights matrix. Row is downstream neurone, column is
    % upstream neurone. Spike raster S (col vect) at a given time step is
    % used to compute PSP change to downstream neurones by N.W * S.
    % Initialise as random value uniformly distributed between 0 and 1.
    % This is the relative strength of each synapse.
    N.W = rand( C.N ) ;
    
    % Multiply relative synapse strength by maximum strength to get final
    % values. Do so for each combination of neurone types.
    N.W( N.e , N.e ) = N.W( N.e , N.e )  .*  C.psp.ee ;
    N.W( N.i , N.e ) = N.W( N.i , N.e )  .*  C.psp.ie ;
    N.W( N.e , N.i ) = N.W( N.e , N.i )  .*  C.psp.ei ;
    N.W( N.i , N.i ) = N.W( N.i , N.i )  .*  C.psp.ii ;
    
    % Neurones are not allowed to synapse onto themselves
    N.W( eye( C.N , 'logical' ) ) = 0 ;
    
    % Store random number state
    N.rng = rng ;
    
    % Return network instance
    varargout = { N } ;
    
    % Restore rng state
    rng( rngtmp )
    
    
  %-- Generate random input currents for excitatory neurones --%
  
  case  'input'
    
    % Set default values for type, par, and rep
    type = 'const' ;
     par =    [ ]  ;
     rep =  false  ;
    
    %- Check input -%
    
    % Check max number of input/output args
     narginchk( 4 , 8 )
    nargoutchk( 0 , 2 )
    
    % Grab first three args, these are now guaranteed
    [ N , off , on ] = varargin{ 1 : 3 } ;
    
    % Basic check on LIF network struct
    if  ~ isstruct( N )
      error( 'lif: input, N must be LIF network struct' )
    end
    
    % Make sure that off and on are scalar, positive values
    if ~isscalar( off ) || ~isnumeric( off ) || ~isfinite( off ) || off < 0
      error( 'lif: input, off must be scalar number >= 0' )
    elseif ~isscalar( on ) || ~isnumeric( on ) || ~isfinite( on ) || on<= 0
      error( 'lif: input, on must be scalar number > 0' )
    end
    
    % Check that there is an even number of name/value arguments, modulus 2
    % returns non-zero value for odd numbers, triggering if statement
    if  mod( nargin - 4 , 2 )
      error( 'lif: input, even number of name/value input arguments' )
    end
    
    % Name of each name/value pair input arguments
    for  arg = 5 : 2 : nargin
      
      % Point to name/value pair
      [ namstr , val ] = varargin{ arg - 1 : arg } ;
      
      % Is namstr even a string?
      if  ~ischar( namstr )  ||  ~isvector( namstr )  ||  ~isrow( namstr )
        error( [ 'lif: input, name of each name/value pair must be ' , ...
          'a string i.e. char row vector' ] )
      end
      
      % Evaluate arguments
      switch  namstr
        
        % Type of current waveform in the ON epoch
        case  { 'const' , 'ramp' , 'sine' , 'noise' , 'norm' }
          
          % Override defaults
          type = namstr ;
           par = val ;
          
          % Check par
          if  ~( isvector( par ) || isempty( par ) ) || ...
               ~isnumeric( par ) || ~all( isfinite( par ) )
            error( [ 'lif: input, par must be a numeric vector of ' , ...
              'finite values, or []' ] )
          end
          
          % Guarantee double floating point
          if  ~ isa( par , 'double' ) , par = double( par ) ; end
          
        % Repeat current time series across LIF neurones?
        case  'repeat'
          
          % Override default
          rep = val ;
          
          % Check rep
          if  ~ isscalar( rep ) || ~ isnumeric( rep ) || ~ isfinite( rep )
            error( 'lif: input, rep must be scalar numeric & finite' )
          end
          
        % Invalid string
        otherwise
          error( 'lif: input, unrecognised name/value pair %s' , namstr )
          
      end % eval args
    end % name/value
    
    % Do not allow rep to be negative when type of signal is not white
    % noise
    if  ~ any( strcmp( type , { 'noise' , 'norm' } ) )  &&  rep < 0
      error( 'lif: input, rep < 0 not defined for type ''%s''' , type )
    end
    
    
    %- Generate input current -%
    
    % Point to network parameters
    C = N.C ;
    
    % Convert from milliseconds to samples, rounding up
    off = ceil( off ./ C.dt ) ;
     on = ceil(  on ./ C.dt ) ;
    
    % Generate off segments
    off = zeros( 1 , off ) ;
    
    % Set random number generator to network's state
    rngtmp = rng( N.rng ) ;
    
    % Select type of current waveform in the ON epoch
    switch  type
      
      % Constant value over time
      case  'const'
        
        % Use default
        if  isempty( par )
          
          par = 1.5 ;
          
        % Otherwise, check that the correct number of values was given
        elseif  numel( par ) ~= 1
          
          error( 'lif: input, const, par must be scalar' )
          
        end % check par
        
        % Current waveform
        I = [ off , par * ones( 1 , on ) , off ] ;
        
      % Ramping current over time
      case   'ramp'
        
        % Use default
        if  isempty( par )
          
          base = 0.0 ;
           amp = 1.5 ;
          
        % Otherwise, check that the correct number of values was given
        elseif  numel( par ) ~= 2
          
          error( 'lif: input, ramp, par must have 2 values' )
          
        % Extract parameters
        else
          
          base = par( 1 ) ;
           amp = par( 2 ) ;
          
        end % check par
        
        % Determine additive constant based on sign of amp. In practice, we
        % only need to add the absolute value of amp onto base if the ramp
        % is decreasing.
        if  amp < 0 , base = base + abs( amp ) ; end
        
        % Build current waveform
        I = [ off , base + ( 1 : on ) ./ on .* amp , off ] ;
        
      % Sinusoidal current oscillation
      case   'sine'
        
        % Use default
        if  isempty( par )
          
          base = 0.0 ;
           amp = 1.5 ;
          freq =  40 ;
          
        % Otherwise, check that the correct number of values was given
        elseif  numel( par ) ~= 3
          
          error( 'lif: input, sine, par must have 3 values' )
          
        % Extract parameters
        else
          
          base = par( 1 ) ;
           amp = par( 2 ) ;
          freq = par( 3 ) ;
          
        end % check par
        
        % Time points, in seconds
        t = ( 1 : on ) .* C.dt ./ 1e3 ;
        
        % Build sinusoidal waveform
        on = base  +  ( 1 + sin( 2*pi * freq * t - pi/2 ) ) ./ 2 .* amp ;
        
        % Build waveform
        I = [ off , on , off ] ;
        
      % Uniformly distributed white noise
      case  'noise'
        
        % Use default
        if  isempty( par )
          
          base = 0.0 ;
           amp = 1.5 ;
          
        % Otherwise, check that the correct number of values was given
        elseif  numel( par ) ~= 2
          
          error( 'lif: input, noise, par must have 2 values' )
          
        % Extract parameters
        else
          
          base = par( 1 ) ;
           amp = par( 2 ) ;
          
        end % check par
        
        % Determine number of unique rows to sample
        if  rep < 0 , rows = C.N ; else , rows = 1 ; end
        
        % Repeat off segment across rows
        off = repmat( off , rows , 1 ) ;
        
        % Generate white noise in ON epoch
        on = ( amp - base ) .* rand( rows , on )  +  base ;
        
        % Build waveform
        I = [ off , on , off ] ;
        
      % Normally distributed white noise
      case   'norm'
        
        % Use default
        if  isempty( par )
          
          avg = 1.5 / 2.0 ;
           sd = avg / 3.5 ;
          
        % Otherwise, check that the correct number of values was given
        elseif  numel( par ) ~= 2
          
          error( 'lif: input, noise, par must have 2 values' )
          
        % Extract parameters
        else
          
          avg = par( 1 ) ;
           sd = par( 2 ) ;
          
        end % check par
        
        % Determine number of unique rows to sample
        if  rep < 0 , rows = C.N ; else , rows = 1 ; end
        
        % Repeat off segment across rows
        off = repmat( off , rows , 1 ) ;
        
        % Generate white noise in ON epoch
        on = sd .* randn( rows , on )  +  avg ;
        
        % Cutoff value above 3.5sd
        cut = avg + 3.5 * sd ;
        
        % Clip values
        on( on > cut ) = cut ;
        
        % Cutoff value below 3.5sd
        cut = avg - 3.5 * sd ;
        
        % Clip values
        on( on < cut ) = cut ;
        
        % Build waveform
        I = [ off , on , off ] ;
            
    end % type of waveform
    
    % Repeat the same time series across neurones
    if  rep > 0 , I = repmat( I , C.N , 1 ) ; end
    
    % Restore rng state and return updated network's state
    N.rng = rng( rngtmp ) ;
    
    %- Done -%
    
    % Return 
    varargout = { I , N } ;
    
    
  %-- Run simulation --%
  
  case  'sim'
    
    % Check max number of input/output args
     narginchk( 3 , 4 )
    nargoutchk( 0 , 2 )
    
    %- Setup -%
    
    % Point to guaranteed input args
    [ N , I ] = varargin{ 1 : 2 } ;
    
    % Basic check on LIF network struct
    if  ~ isstruct( N )
      error( 'lif: sim, N must be LIF network struct' )
    end
    
    % Point to network parameters
    C = N.C ;
    
    % Check input current format
    if  ~ isa( I , 'double' )
      
      error( 'lif: sim, I must be double floating point' )
      
    % No input
    elseif  isempty( I )
      
      error( 'lif: sim, I is empty' )
      
    % This is a row vector or matrix
    elseif  ~ isrow( I )  &&  ~ ismatrix( I )
      
      error( 'lif: sim, I must be a row vector or matrix' )
      
    % Check number of rows match number of neurones, if matrix
    elseif  ~ isvector( I )  &&  ismatrix( I )  &&  C.N ~= size( I , 1 )
      
      error( 'lif: sim, I must have one row per LIF neurone' )
      
    % Invalid numerical values
    elseif  ~ all( isfinite( I ) , 'all' )
      
      error( 'lif: sim, I must have finite numerical values' )
      
    end % check I
    
    % Has volflg been provided?
    if  nargin == 4
      
      % Get it
      volflg = varargin{ 3 } ;
      
      % Check volflg
      if  ~ isscalar( volflg )  ||  ~( isnumeric( volflg ) || ...
          islogical( volflg ) )
        
        error( 'lif: sim, volflg must be scalar logical or numeric' )
        
      end % check
      
    % No, set default
    else
      
      volflg = true ;
      
    end % volflg
    
    % Convert unit of input from current in nA to voltage in mV. Now the I
    % simply stands for 'I'nput.
    I = I ./ C.C ;
    
    % Number of time steps
    Nt = size( I , 2 ) ;
    
    % Get input current scaling value for each neurone, by type
    scale = zeros( C.N , 1 ) ;
    scale( N.e ) = C.Ie ;
    scale( N.i ) = C.Ii ;
    
    % Allocate output
    S.spk = false( C.N , Nt ) ;
    
    if  volflg
      S.vol = zeros( C.N , Nt ) ;
    else
      S.vol = [ ] ;
    end
    
    % Store state of random number generator and apply network's state
    rngtmp = rng( N.rng ) ;
    
    % Randomly initialise membrane voltage of all neurones on first time
    % step, uniformly distributed from resting voltage to spiking
    % threshold.
    V = ( C.V.threshold - C.V.leak ) .* rand( C.N , 1 )  +  C.V.leak ;
    
    % Add input on time step 1
    V = V  +  scale .* I( : , 1 ) ;
    
    % Find spiking units, this caries forward to hyperpolarisation step
    i = V  >=  C.V.threshold ;
    
    % Set their voltage to spiking level
    V( i ) = C.V.spike ;
    
    % Store initial state
    S.spk( i , 1 ) = 1 ;
    if  volflg , S.vol( : , 1 ) = V ; end
    
    % And also store the current state of the random number generator
    N.rng = rng ;
    
    % Restore rng state
    rng( rngtmp )
    
    
    %- Run simulation -%
    
    % Time, excluding first step
    for  t = 2 : Nt
      
      % Apply membrane current leak
      V = V  +  C.dt .* ( -( V - C.V.leak ) ./ C.tau ) ;
      
      % Add external input current to neurones
      V = V  +  scale .* I( : , t ) ;
      
      % Hyperpolarise neurones that fired at t - 1
      V( i ) = C.V.reset ;
      
      % Add post-synaptic potentials from incoming spikes
      V = V  +  N.W * i ;
      
      % Find spikes that fire on current time step, used in next iteration
      i = V  >=  C.V.threshold ;
      
      % Set spiking neurones' membrane voltage to spiking level
      V( i ) = C.V.spike ;
      
      % Store current state of network
      S.spk( i , t ) = 1 ;
      if  volflg , S.vol( : , t ) = V ; end
      
    end % time
    
    %- Done -%
    
    % Return simulation results
    varargout = { S , N } ;
    
    
  %-- Function string is not recognised --%
  
  otherwise , error( 'fstr string unrecognised: %s' , fstr )
  
    
end % func selection

