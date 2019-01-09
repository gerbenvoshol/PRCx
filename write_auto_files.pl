#!/usr/bin/perl
#   
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#                  PRC, the profile comparer, version 1.5.6
#
#   write_auto_files.pl: Perl script for generating the auto_DP_xxx.c files
#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#  Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
#  Copyright (C) 2018-19 Gerben Voshol, Leiden University, The Netherlands
#  All Rights Reserved
#
#  This source code is distributed under the terms of the GNU General Public 
#  License. See the files COPYING and LICENSE for details.
#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#
#  This script automatically generates the following files:
#
#     auto_DP_comps.c
#     auto_DP_inner_loop.c
#     forward_transitions.txt
#     backward_transitions.txt
#
#  General conventions and typical values:
#
#    'M' =  match state
#    'I' = insert state
#    'D' = delete state
#
#    'x' = no transition
#    'o' = transition from the same node
#    '-' = transition from the previous node
#    '+' = transition from the following node
#
#    $dst     = 'MD'
#    $dst1    = 'M'
#    $diff    = 'x-'
#    $diff1   = 'x'   
#    $src     = 'MM'
#    $src1    = 'M'
#    $source  = 'MxM-'
#    $source1 = 'Mx' 
#
#    %prof_trans = hash of allowed transitions, as below
#    $prof_trans = 'PLAN7'
#
use strict;

our (%algos, %prof_trans, %pair_states, %pair_trans, %want_comps);

%algos = (  'VITERBI' =>   {  'DIRECTION'  => 'DIR_FORWARD', 
			      'ALGO_PATHS' => 'BEST_PATH'  },
	    
	    'VIT_BACK' =>  {  'DIRECTION'  => 'DIR_BACKWARD',
			      'ALGO_PATHS' => 'BEST_PATH'  },
	    
	    'FORWARD' =>   {  'DIRECTION'  => 'DIR_FORWARD',
			      'ALGO_PATHS' => 'SUM_PATHS'  },
	    
	    'BACKWARD' =>  {  'DIRECTION'  => 'DIR_BACKWARD',
			      'ALGO_PATHS' => 'SUM_PATHS'  },
	    
	    'OPT_ACCUR' => {  'DIRECTION'  => 'DIR_FORWARD',
			      'ALGO_PATHS' => 'BEST_PATH'  }  );

#                       dst     .................. sources .................
$prof_trans{PLAN7} = {  'M' => { 'Mx' => 1, 'M-' => 1, 'I-' => 1, 'D-' => 1 },
			'I' => { 'Ix' => 1, 'Mo' => 1, 'Io' => 1            },
			'D' => { 'Dx' => 1, 'M-' => 1,            'D-' => 1 }
		     };

$prof_trans{PLAN9} = {  'M' => { 'Mx' => 1, 'M-' => 1, 'I-' => 1, 'D-' => 1 },
			'I' => { 'Ix' => 1, 'Mo' => 1, 'Io' => 1, 'Do' => 1 },
			'D' => { 'Dx' => 1, 'M-' => 1, 'I-' => 1, 'D-' => 1 }
		     };

$pair_states{SPACE5} 
= { 'MM'=>1, 'MI'=>1, 'IM'=>1, 'MD'=>1, 'DM'=>1 };

$pair_states{SPACE6}
= { 'MM'=>1, 'MI'=>1, 'IM'=>1, 'MD'=>1, 'DM'=>1, 'DD'=>1 };

$pair_states{SPACE8}
= { 'MM'=>1, 'MI'=>1, 'IM'=>1, 'MD'=>1, 'DM'=>1, 'DD'=>1, 'DI'=>1, 'ID'=>1 };

$pair_states{SPACE9}
= { 'MM'=>1, 'MI'=>1, 'IM'=>1, 'II'=>1, 
    'MD'=>1, 'DM'=>1, 'DD'=>1, 'DI'=>1, 'ID'=>1 };

&populate_pair_trans;
&print_comps;
&print_search;
&print_which_comps;
exit;


# initialize %pair_trans with allowed transitions in the pair HMM
#
sub populate_pair_trans
{
    my ($prof_trans, $pair_states, $pair_trans);

    # keep track of approved pair HMM transitions
    open F, '>forward_transitions.txt' or die;
    open B, '>backward_transitions.txt' or die;
    print F &banner;
    print B &banner;

    for $prof_trans  ('PLAN7', 'PLAN9')                        {
    for $pair_states ('SPACE5', 'SPACE6', 'SPACE8', 'SPACE9' ) {
    for $pair_trans  ('VIA_MM', 'ALL_TRANS')
    {
	print F "$prof_trans $pair_states $pair_trans:\n";
	print B "$prof_trans $pair_states $pair_trans:\n";
	my $total = 0;

	for my $dst   ( sort keys %{ $pair_states{$pair_states} } ) {    
	for my $src   ( sort keys %{ $pair_states{$pair_states} } ) {
	for my $diff1 ( 'x', 'o', '-' ) {
	for my $diff2 ( 'x', 'o', '-' )
	{
	    # bureaucracy
	    my ($dst1,  $dst2)  = split //, $dst;
	    my ($src1,  $src2)  = split //, $src;
	    my $diff    = "$diff1$diff2";
	    my $source1 = "$src1$diff1";
	    my $source2 = "$src2$diff2";
	    my $source  = "$source1$source2";
	    
	    # if 'VIA_MM', all transitions (apart from transitions to self) 
	    # must go via MM
	    next if $pair_trans eq 'VIA_MM' and
		$src ne 'MM' and $dst ne 'MM' and $src ne $dst;
	    
	    # check that the transitions are OK with the profile HMMs
	    next unless exists $prof_trans{$prof_trans}{$dst1}{$source1};
	    next unless exists $prof_trans{$prof_trans}{$dst2}{$source2};

	    # check that the diffs are sound
	    if   ( $dst eq 'DD' ) { next unless $diff eq '--' }
	    elsif( $dst1 eq 'D' ) { next unless $diff eq '-x' }
	    elsif( $dst2 eq 'D' ) { next unless $diff eq 'x-' }
	    else                  { next if $diff1 eq 'x' or $diff2 eq 'x' };
	    
	    # IoIo -> II is horrible and nasty
	    next if $dst eq 'II' and $source eq 'IoIo';

	    # all clear, transition approved
	    $pair_trans{$prof_trans}{$pair_states}{$pair_trans}
	    {'DIR_FORWARD'}{$dst}{$source} = 1;
	    print F "$source -> $dst\n";

	    # so far this was a DIR_FORWARD transition
	    # now calculate the corresponding DIR_BACKWARD transition
	    my $back_dst    = $src;
	    my $back_diff1  = $diff1 eq '-' ? '+' : $diff1;
	    my $back_diff2  = $diff2 eq '-' ? '+' : $diff2;
	    my $back_source = "$dst1$back_diff1$dst2$back_diff2";

	    $pair_trans{$prof_trans}{$pair_states}{$pair_trans}
	    {'DIR_BACKWARD'}{$back_dst}{$back_source} = 1;
	    print B "$back_source -> $back_dst\n";

	    $total++; 
	}}}};

	print F "(total = $total)\n\n";
	print B "(total = $total)\n\n";
    }}};

    close F;
    close B;
};


# write the definitions of compX_sum and compX_best to auto_DP_comps.c
#
sub print_comps
{
    my (@lines);


    push @lines, '#if ALGO_SPACE==LINSPACE', '';

    push @lines,
    (  
       '#ifdef WANT_COMP1_SUM',
       '#error "VIA_MM will not work with SPACE8 or SPACE9!"',
       '#endif',
       ''
    );

    for my $n (2 .. 9)
    {
	push @lines,
	(  
	   "#ifdef WANT_COMP$n\_SUM",
	   "#define comp$n\_sum comp$n\_sum\_lin",
	   "static double comp$n\_sum_lin ( "
	   .  join( ', ', map("double sc$_", 1 .. $n) )  .  ' )',
	   '{',
	   'return '  .  join( '+', map("sc$_", 1 .. $n) )  .  ';',
	   '}',
	   '#endif',
	   ''
	);
    };

    push @lines, '#endif', '', '';
    push @lines, '#if ALGO_SPACE==LOGSPACE', '';
    
    push @lines, 
    ( 
       '#ifdef WANT_COMP1_SUM',
       '#error "VIA_MM will not work with SPACE8 or SPACE9!"',
       '#endif',
       '',
       "#ifdef WANT_COMP2_SUM",
       '#define comp2_sum comp2_sum_log',
       'static double comp2_sum_log ( double sc1, double sc2 )',
       '{',
       'double sc;',
       '',
       'if( sc1 < sc2 ) sc = sc2 + log(1.0 + exp(sc1-sc2));',
       'else            sc = sc1 + log(1.0 + exp(sc2-sc1));',
       '',        
       'return sc;',
       '}',
       '#endif',
       ''
    );

    for my $n (3 .. 9)
    {
	push @lines,
	(  
	   "#ifdef WANT_COMP$n\_SUM",
	   "#define comp$n\_sum comp$n\_sum\_log",
	   "static double comp$n\_sum_log ( "
	   .  join( ', ', map("double sc$_", 1 .. $n) )  .  ' )',
	   '{',
	   'double max = sc1;',
	   '',
	   map("if( max < sc$_ ) max = sc$_;", 2 .. $n), 
	   '',
	   'return max+log('  
	   .  join( '+', map("exp(sc$_-max)", 1 .. $n) )  .  ');',
	   '}',
	   '#endif',
	   ''
	);
    };

    push @lines, '#endif', '', '';
    

    push @lines,
    (  
       '#ifdef WANT_COMP1_BEST',
       '#error "VIA_MM will not work with SPACE8 or SPACE9!"',
       '#endif'
    );

    for my $n (2 .. 9)
    {
	push @lines,
	(  
	   "#ifdef WANT_COMP$n\_BEST",

	   # header
	   "static double comp$n\_best ( int *tr, "
	   .  join( ', ', map("int st$_, double sc$_", 1 .. $n) )  .  ' )',

	   # body
	   '{',
	   'double sc = sc1;',
	   'int   st = st1;',
	   '',
	   map("if( sc$_ > sc ){ sc = sc$_; st = st$_; };", 2 .. $n),
	   '',
	   '*tr = st;',
	   'return sc;',
	   '}',

	   '#endif',
	   ''
        );
    };

    open F, '>auto_DP_comps.c' or die;
    print F &indent(0, @lines);
    close F;
};

# write the inner loops to auto_DP_inner_loop.c
# 
sub print_search
{
    my (@lines);
     
    for my $prof_trans  ('PLAN7', 'PLAN9')                              {
    for my $pair_states ('SPACE5', 'SPACE6', 'SPACE8', 'SPACE9')        {
    for my $pair_trans  ('VIA_MM', 'ALL_TRANS')                         {
    for my $algorithm   ('VITERBI', 'VIT_BACK', 
			 'FORWARD', 'BACKWARD', 'OPT_ACCUR') 
    {
	# params = topology & algorithm
	my $params = { prof_trans  => $prof_trans,
		       pair_states => $pair_states,
		       pair_trans  => $pair_trans,
		       algorithm   => $algorithm };

	# get the order in which the states will be calculated
	my @dsts;

	if( $pair_states eq 'SPACE9' )
	{
	    # for SPACE9 need to be careful about II
	    if( $algorithm eq 'VITERBI' or $algorithm eq 'FORWARD' )
	    {
		# II depends on the other states
		@dsts = ( sort(keys %{ $pair_states{'SPACE8'} }), 'II' );
	    }
	    elsif( $algorithm eq 'OPT_ACCUR' )
	    {
		# no II for OPT_ACCUR!
		@dsts = sort(keys %{ $pair_states{'SPACE8'} });
	    }
	    elsif( $algorithm eq 'VIT_BACK' or $algorithm eq 'BACKWARD' )
	    {
		# other states depend on II
		@dsts = ( 'II', sort(keys %{ $pair_states{'SPACE8'} }) );
	    };
	}
	else
	{
	    # SPACE8 and less - the order is irrelevant
	    @dsts = sort keys %{ $pair_states{$pair_states} };
	};


	push @lines, &params_if_endif
	    ( $params, map &dst_recursion($params, $_), @dsts );

    }}}};
    
    open F, '>auto_DP_inner_loop.c' or die;
    print F &indent(1, @lines);
    close F;
};

# write the WANT_COMPX_XXXX requirements into auto_DP_which_comps.c
# 
sub print_which_comps
{
    my (@lines);
     
    for my $prof_trans  ('PLAN7', 'PLAN9')                              {
    for my $pair_states ('SPACE5', 'SPACE6', 'SPACE8', 'SPACE9')        {
    for my $pair_trans  ('VIA_MM', 'ALL_TRANS')                         {
    for my $algorithm   ('VITERBI', 'VIT_BACK', 
			 'FORWARD', 'BACKWARD', 'OPT_ACCUR') 
    {
	# params = topology & algorithm
	my $params = { prof_trans  => $prof_trans,
		       pair_states => $pair_states,
		       pair_trans  => $pair_trans,
		       algorithm   => $algorithm };

	push @lines, &params_if_endif
	    (  $params,
	       map ( '#define WANT_'.uc($_), 
		     sort keys %{ $want_comps
				  {$prof_trans}
				  {$pair_states}
				  {$pair_trans}
				  {$algorithm} } ),
	       '' );

    }}}};
    
    open F, '>auto_DP_which_comps.c' or die;
    print F &indent(0, @lines);
    close F;
};

# return something like
#
#   (  "S[m1][m2][sMI] = X[m1][m2][xMI] x",
#      "comp2_sum",
#      "(",
#      "S[m1-1][m2][sMI] x St1[m1][tMM] x St2[m2][tII],",
#      "S[m1-1][m2][sMM] x St1[m1][tMM] x St2[m2][tMI]",
#      ");"
#   )
# 
# i.e. the full recursion for S[m1][m2][s$dst] 
#
sub dst_recursion
{
    #            $dst is the short form (MD, DI etc.)
    my ($params, $dst) = @_;

    my @sources = sort keys %{ $pair_trans
			       {$params->{prof_trans}}
			       {$params->{pair_states}}
			       {$params->{pair_trans}}
			       {$algos{$params->{algorithm}}{'DIRECTION'}}
			       {$dst} };

    # ignore II for OPT_ACCUR
    if( $params->{algorithm} eq 'OPT_ACCUR' )
    {
	for my $i (0 .. $#sources)
	{
	    splice(@sources, $i, 1) if $sources[$i] =~ /^I.I.$/;
	};
    };
    
    my $algorithm  = $params->{'algorithm'};
    my $algo_paths = $algos{$algorithm}{'ALGO_PATHS'};
    my $n          = scalar(@sources);

    my ($comp_func, @fn, @lines);


    if( $algo_paths eq 'BEST_PATH' )
    {
	$comp_func = "comp$n\_best";
	&want_comp($params, $comp_func);
	push @fn, $comp_func, '(', "&Tr[m1][m2][s$dst],";

	for my $source (@sources)
	{
	    my $rec = &source_recursion($algorithm, $dst, $source);
	    my ($src1, $diff1, $src2, $diff2) = split //, $source;
	    push @fn, "s$src1$src2, $rec,";
	};
    }
    elsif( $algo_paths eq 'SUM_PATHS' )
    {
	$comp_func = "comp$n\_sum";
	&want_comp($params, $comp_func);
	push @fn, $comp_func, '(';

	for my $source (@sources)
	{
	    my $rec = &source_recursion($algorithm, $dst, $source);
	    push @fn, "$rec,";
	};
    };

    # cut the comma after the last parameter
    $fn[-1] =~ s/\,$//;


    if( $algorithm eq 'OPT_ACCUR' )
    {
	&want_comp($params, 'comp2_sum');

	push @lines, 
	(  
	   "S[m1][m2][s$dst] = ",
	   'comp2_sum',
	   '(',
	   "S[m1][m2][s$dst],",
	   @fn,
	   ')',
	   ');',
	   ''
	);

	push @lines, 
	(  
	   '#if ALGO_SPACE == LINSPACE',
	   #'if( (S[m1][m2][sMM] > dp->p->max_double) || '
	   #. '( (S[m1][m2][sMM]!=0.0) && '
	   #. '(S[m1][m2][sMM] < dp->p->min_double) ) )',
	   #'dp->lin_overflow=1;',
	   'if( S[m1][m2][sMM] > MAX_FLOAT ) ',
	   '{',
	   'dp->lin_overflow=1;',
	   '};',
	   '#endif',
	   '',
	   'if( dp->best_S < (temp = S[m1][m2][sMM] x StOut) )',
	   '{',
	   'dp->best_S  = temp;',
	   'dp->best_m1 = m1;',
	   'dp->best_m2 = m2;',
	   '};',
	   '',
	)
	    if $dst eq 'MM';
    }
    elsif( $algo_paths eq 'BEST_PATH' )
    {
	# VITERBI or VIT_BACK

	if( $dst eq 'MM' )
	{
	    push @lines,
	    (  
	       'S[m1][m2][sMM] = ', 
	       @fn, 
	       ');', 
	       '',
	       '#if ALGO_SPACE == LINSPACE',
	       #'if( (S[m1][m2][sMM] > dp->p->max_double) || '
	       #. '( (S[m1][m2][sMM]!=0.0) && '
	       #. '(S[m1][m2][sMM] < dp->p->min_double) ) )',
	       #'dp->lin_overflow=1;',
	       'if( S[m1][m2][sMM] > MAX_FLOAT ) ',
	       '{',
	       'dp->lin_overflow=1;',
	       '};',
	       '#endif',
	       '',
	       'if( dp->best_S < (temp = S[m1][m2][sMM] x StOut) )',
	       '{',
	       'dp->best_S  = temp;',
	       'dp->best_m1 = m1;',
	       'dp->best_m2 = m2;',
	       '};',
	       '',
	       'S[m1][m2][sMM] xeq X[m1][m2][xMM];',
	       '',
	       'if( S[m1][m2][sMM] < StIn )',
	       '{',
	       'S[m1][m2][sMM] = StIn;',
	       '',
	       '// N.B. sBE is now indicated by Tr < 0',
	       'Tr[m1][m2][sMM] -= 20;',
	       '};',
	       ''
	    );
	}
	elsif( $dst eq 'MI' or $dst eq 'IM' or $dst eq 'II' )
	{
	    my ($dst1, $dst2) = split //, $dst;

	    push @lines, 
	    (  
	       "S[m1][m2][s$dst] = X[m1][m2][x$dst] x ", 
	       @fn,
	       ');',
	       ''
	    );
	}
	else
	{
	    push @lines, 
	    (  
	       "S[m1][m2][s$dst] = ", 
	       @fn, 
	       ');', 
	       ''
	    );
	};
    }
    else
    {
	# FORWARD or BACKWARD

	if( $dst eq 'MM' )
	{
	    &want_comp($params, 'comp2_sum');

	    push @lines, 
	    (  
	       'S[m1][m2][sMM] = ',
	       @fn, 
	       ');',
	       '',
	       'other_S[m1][m2][sMM] xeq S[m1][m2][sMM];',
	       '',
	       '#if ALGO_SPACE == LINSPACE',
	       #'if( (S[m1][m2][sMM] > dp->p->max_double) || '
	       #. '( (S[m1][m2][sMM]!=0.0) && '
	       #. '(S[m1][m2][sMM] < dp->p->min_double) ) )',
	       #'dp->lin_overflow=1;',
	       'if( S[m1][m2][sMM] > MAX_FLOAT ) ',
	       '{',
	       'dp->lin_overflow=1;',
	       '};',
	       '#endif',
	       '',
	       'temp = S[m1][m2][sMM] x StOut;',
	       'dp->sum_S = comp2_sum( dp->sum_S, temp );',
	       '',
	       'if( dp->best_S < temp )',
	       '{',
	       'dp->best_S  = temp;',
	       'dp->best_m1 = m1;',
	       'dp->best_m2 = m2;',
	       '};',
	       '',
	       'S[m1][m2][sMM] = ',
	       'comp2_sum', 
	       '(',
	       'S[m1][m2][sMM] x X[m1][m2][xMM],',
	       'StIn',
	       ');'
	    );
	}
	elsif( $dst eq 'MI' or $dst eq 'IM' or $dst eq 'II' )
	{
	    my ($dst1, $dst2) = split //, $dst;

	    push @lines,
	    (  
	       "S[m1][m2][s$dst] = ", 
	       @fn,
	       ');',
	       '',
	       "other_S[m1][m2][s$dst] xeq S[m1][m2][s$dst];",
	       '',
	       "S[m1][m2][s$dst] xeq X[m1][m2][x$dst];",
	       ''
	    );
	}
	else
	{
	    push @lines, 
	    (  
	       "S[m1][m2][s$dst] = ", 
	       @fn, 
	       ');', 
	       '',
	       "other_S[m1][m2][s$dst] xeq S[m1][m2][s$dst];",
	       ''
	    );
	};
    };

    return @lines;
};


# return something like
#
#   "S[m1-1][m2][sDM] x St1[m1][tDD]"
#
sub source_recursion
{
    my (  # 'VITERBI' or 'OPT_ACCUR' etc.
	  $algorithm,

	  # MI, DD etc.
	  $dst,
 
	  # for forward transitions:  MxD-, M-Mo etc.
	  # for backward transitions: M+M+, MxD+ etc. 
	  $source

        ) = @_;

    my ($src1, $diff1, $src2, $diff2) = split //, $source;
    my ($dst1, $dst2)                 = split //, $dst;

    my ($src_ind1, $src_ind2, $tr_ind1, $tr_ind2, $line);

    if( $algos{$algorithm}{'DIRECTION'} eq 'DIR_FORWARD' )
    {
	$src_ind1 = $diff1 eq '-' ? 'm1-1' : ' m1 ';
	$src_ind2 = $diff2 eq '-' ? 'm2-1' : ' m2 ';
	$tr_ind1  = 'm1';
 	$tr_ind2  = 'm2';
    }
    elsif( $algos{$algorithm}{'DIRECTION'} eq 'DIR_BACKWARD' )
    {
	$src_ind1 = $diff1 eq '+' ? 'm1+1' : ' m1 ';
	$src_ind2 = $diff2 eq '+' ? 'm2+1' : ' m2 ';
	$tr_ind1  = $src_ind1;
	$tr_ind2  = $src_ind2;
    };

    $line  = "S[$src_ind1][$src_ind2][s$src1$src2]";

    if( $algorithm ne 'OPT_ACCUR' )
    {
	if( $algos{$algorithm}{'DIRECTION'} eq 'DIR_FORWARD' )
	{
	    $line .= " x St1[$tr_ind1][t$src1$dst1]" unless $diff1 eq 'x';
	    $line .= " x St2[$tr_ind2][t$src2$dst2]" unless $diff2 eq 'x';    
	}
	elsif( $algos{$algorithm}{'DIRECTION'} eq 'DIR_BACKWARD' )
	{
	    $line .= " x St1[$tr_ind1][t$dst1$src1]" unless $diff1 eq 'x';
	    $line .= " x St2[$tr_ind2][t$dst2$src2]" unless $diff2 eq 'x';
	};
    };

    return $line;
};

# record that we want $comp_func for $params->{algorithm} under a given
# topology 
#
sub want_comp
{
    my ($params, $comp_func) = @_;

    $want_comps
    {$params->{prof_trans}}
    {$params->{pair_states}}
    {$params->{pair_trans}}
    {$params->{algorithm}}
    {$comp_func} = 1;
};

# wrap @lines in #if ... #endif based on $params
#
sub params_if_endif
{
    my $params = shift @_;

    return ( "#if PROF_HMM_TRANS  == $params->{prof_trans}",
	     "#if PAIR_HMM_STATES == $params->{pair_states}",
	     "#if PAIR_HMM_TRANS  == $params->{pair_trans}",
	     "#if ALGORITHM       == $params->{algorithm}",
	     '',
	     @_,
	     '#endif',
	     '#endif',
	     '#endif',
	     '#endif',
	     '', 
	     '' );	     
};

sub indent
{
    my ($n_tabs, @code) = @_;
    my ($no_adjust, @lines);

    foreach my $line (@code)
    {
	unless($line =~ /^\#/)
	{
	    # should the next line be indented more, the same, or less?
	    my $diff = ($line =~ s/([({])/$1/g) - ($line =~ s/([)}])/$1/g);

	    # ) and } need to be indented the same way as the next line
	    $n_tabs += $diff if $diff < 0;

	    # note the $no_adjust part (for lines that follow a backslash)
	    $line = ("\t" x $n_tabs) . $line unless $no_adjust;
	    
	    # finish adding the diff
	    $n_tabs += $diff if $diff > 0;
	};

	if( $line =~ s/\\$// ){ $no_adjust = 1 }
	else   { $line .= "\n"; $no_adjust = 0 };

	push @lines, $line;
    };

    return ( &banner, @lines );
};


sub banner
{
    return "// automatically generated by write_DP_loops.pl\n\n";
};
