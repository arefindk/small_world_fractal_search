: # Find perl...
eval 'exec perl5 -w -S $0 "$@"' 
if 0;

#!/usr/bin/perl5 -w
#
############################################################
# $Id: DiskFracDim.pm,v 1.4 2000/10/24 18:00:46 lw2j Exp $
# 
# $Log:	DiskFracDim.pm,v $
# Revision 1.4  2000/10/24  18:00:46  lw2j
# Now supports $q=1.  Logarithms now base 2.
# 
# Revision 1.3  2000/10/06  19:19:44  lw2j
# Fixed bugs involving q != 2.  Barred q=1, for now.
# 
# Revision 1.2  2000/10/03  14:08:01  lw2j
# Added magic perl invoker.
# 
# Revision 1.1  2000/07/18  13:27:11  lw2j
# Initial revision
#   
############################################################
#
# Purpose:  
#  To provide simple code for computing the 'correlation' fractal
#  dimension $D_{q}$.
#
#  It's {\em very} simple code, so the data is expected to be accessed 
#  via a wrapper, for reasons of flexibility.  This code is 
#
#                   ***** NOT *****
#
#  optimized for speed in an industrial sense.
#
#  That is, the module accepts as input an object reference.  This
#  object reference must support:
#
#     get_count($)       # return how many items there are
#     get_dims($)        # return embedding dimensionality
#     get_obj($$)        # return the object by 0-based index, as an
#                        # array of scalars
# 
#  The object should not change while we're processing it...
#
#  This object ref will be 'data_o' below.
# 
#  Other parameters that should be worried about:
#
#     method             # 'fast' (inexact box-count) or 'slow'
#                        # (iterate over pairwise distances, but
#                        #  more exact; DO NOT USE if you care
#                        # about speed)
#     r_min              # minimum radius for our approximation
#     r_max              # maximum radius for our approximation
#     r_count            # maximum number of radii to use (>= 1)
#     intervals          # if defined (default 20) and the method 
#                        # is 'slow', the log-distance axis is
#                        # divided into this many intervals, and
#                        # the point with the LOWEST pair-count is
#                        # chosen for the line-fitting if any 
#                        # exists.
#     q                  # normally 2, the degree
#
#
#     verbose            # if nonzero, log radius/log count values
#                        # are printed, prepended by ';;; '.
#
#  The default values tend to be fairly sane, but if the data 
#  you have has ranges significantly different than mine, you
#  may want to adjust the radii (for instance).  See 
#  'set_default_params'.
#
#  Logarithms, incidentally, should now be base 2 throughout.
#  In cases where the scaling applies to both axes, it does 
#  not matter from the point of determining fractal dimension --
#  since it does not affect the slope -- but it does matter
#  in terms of absolute numbers, and for $q=1 (based on 
#  entropy, which normally uses bits).
#
#  The methods of greatest concern are new, set_params, clear
#  and fracdim.
#  
#  Sample usage:   
#    See 'FracTest.pl', which attempts to load a data file 'tri'.
#
##########################################################################
#
# Package Header
#
##########################################################################

package DiskFracDim;
require Exporter;


##### <change_me>
#
# Change this to where you placed the modules.
#
use lib '/afs/cs.cmu.edu/user/christos/TST/FracDim/FracDim';
##### </change_me>

use strict;
use DB_File;
use POSIX;
use mktemp;      # temporary file support (for boxcount database)

##########################################################################

BEGIN {
    use Exporter ();
    use Symbol;
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

    $VERSION     = 1.0;
    @ISA         = qw(Exporter);
    @EXPORT =  qw(
                  new
                  clear
                  set_params
                  fracdim

		  set_default_params
                  
                  compute_ranges
		  box_put
                  box_count
                  all_boxes
		  all_boxes_int
		  all_pairs

		  print_loglog

		  filter_points
                  trim_ranges
                  lsfit
                 );
    %EXPORT_TAGS = ();
    @EXPORT_OK   = ();

    return 1;
}



# new:  creates blessed FracDim objects with default params
#
# Standard usage:
#   
# my $fracObj = FracDim::new();
#
# Extended usage:
#
# my $fracObj = FracDim::new( param => value, param => value, ... ).
sub new($@) {
    my $class = shift;
    my $self  = +{};

    bless $self, $class;

    $self->clear();
    $self->set_default_params();
    $self->set_params(@_);

    # An on-disk btree is used to store box-counts in the 'fast'
    # method -- it has to be on-disk, because in higher 
    # dimensionalities with large data sets the size of the hash
    # table could get *very* large.  

    $DB_BTREE->{cachesize}=100000;

    return $self;   
}



sub set_default_params($) {
    my $self = shift;

    # default:  2^-20 to 2^18, should result in a multiplier of 2
    #
    # I *highly* advise choosing values that give integral 
    # multipliers, as various optimizations can't happen
    # otherwise.
    $self->{'params'} = 
        +{
            'q'       => 2,        
	    normalize => 0,
            r_min     => 0.00000095367431640625,
            r_max     => 262144,
            r_count   => 39,     
	    intervals => 20,
	    verbose   => 0,
	    method    => 'fast',
	    data_o    => undef   
         };
}


# Forget about old data.  If you only call 'new' once, and 
# use set_params to switch data objects, you'll want to
# clear() the old data out of there.
#
# Standard usage:
# 
# $fdObj->clear();
sub clear($) {
    my $self = shift;

    # Forget about the old data.  
    if (exists($self->{'data'})) {
	if (exists($self->{'data'}->{'old_dbref'})) {
	    delete $self->{'data'}->{'old_dbref'};
	    untie %{$self->{'data'}->{'old_boxdata'}};
	    untemp $self->{'data'}->{'old_tnam'};
	}
    }

    $self->{'data'}       = +{};
   
    return $self;
}




# Set zero or more parameters, doing some minimal sanity checking.
#
# Standard usage:
#
# $fdObj->set_params('q' => 2, r_count => 200);
sub set_params($@) {
    my $self  = shift;
    my $param = undef;
    my $arg   = undef;

    while (defined($param = shift)) {
	$arg = shift;

	if (!exists($self->{'params'}->{$param})) {
            die "Unknown param/arg pair '$param':'$arg'";
        }

        $self->{'params'}->{$param} = $arg;
    }

    if (exists($self->{'params'}->{'method'})) {
	my $method = $self->{'params'}->{'method'};

	if (($method ne 'fast') && ($method ne 'slow')) {
	    die "Parameter 'method' must be either 'fast' or 'slow'.";
	}

	if ($method eq 'slow') {
	    if (defined($self->{'params'}->{'intervals'}) &&
		($self->{'params'}->{'intervals'} < 10)) {
		
		# Lower this threshold if you REALLY want to.  Do NOT use
		# just 1, however.
		die "Using below 10 intervals is REALLY not recommended if you're doing the full pairwise method anyway.";
	    }
	}
    } else {
	die "The 'method' parameter is mandatory.";
    }

    return $self;
}


# Returns the $D_{q}$ fractal dimension, the y-intercept, and
# the correlation coefficient.  It's really just a macro.
#
# The data is not automatically clear()-ed afterwards.
sub fracdim($) {
  my $self  = shift;
  my $y_int = undef;
  my $slope = undef;
  my $corr  = undef;

  $self->compute_ranges();

  if ($self->{'params'}->{'method'} eq 'fast') {
      $self->all_boxes();
  } elsif ($self->{'params'}->{'method'} eq 'slow') {
      $self->all_pairs();
  } else {
      die "You really shouldn't reach this point.";
  }

  $self->trim_ranges();

  if ($self->{'params'}->{'method'} eq 'slow') {
      if (defined($self->{'params'}->{'intervals'})) {
	  $self->filter_points();
      }
  }


  if ($self->{'params'}->{'verbose'}) {
      $self->print_loglog();
  }

  ($slope, $y_int, $corr) = $self->lts();

  my $q = $self->{'params'}->{'q'};

  if ($q != 1) {
      $slope /= ($q-1);
  }

  return ($slope, $y_int, $corr);
}


# Compute the minimum and maximum for each dimension. 
# It's required for normalization purposes.
#
# Standard usage:
# 
# $fdObj->compute_ranges();
sub compute_ranges($) {
    my $self   = shift;
    my $dims   = undef;
    my $data_o = undef;
    my $count  = undef;
    my $obj_idx = undef;
    
    my @lows    = ();
    my @highs   = ();

    if (!exists($self->{'params'}->{'data_o'})) {
        die "Corrupted object, doesn't even have an undef data_o";
    }

    if (!exists($self->{'params'}->{'data_o'})) {
        die "Incomplete object, needs a defined data_o reference.";
    }

    $data_o = $self->{'params'}->{'data_o'};
    $dims   = $data_o->get_dims();
    $count  = $data_o->get_count();

    $self->{'data'}->{'dims'}  = $dims;
    $self->{'data'}->{'count'} = $count;

    for ($obj_idx = 0; $obj_idx < $count; $obj_idx++) {
        my @obj_data = undef;
        my $dim_idx  = undef;

        @obj_data = $data_o->get_obj($obj_idx);

        defined(@obj_data) || die "Failed to retrieve data obj $obj_idx";
        (scalar(@obj_data) == $dims) || die "Unexpected size for obj $obj_idx";

	if ($obj_idx == 0) {
	    @lows = @obj_data;
	    @highs = @obj_data;
	} else {
	    for ($dim_idx=0; $dim_idx < $dims; $dim_idx++) {
		if ($lows[$dim_idx] > $obj_data[$dim_idx]) {
                    $lows[$dim_idx] = $obj_data[$dim_idx];
                }
		
                if ($highs[$dim_idx] < $obj_data[$dim_idx]) {
                    $highs[$dim_idx] = $obj_data[$dim_idx];
                }
            }
        }
    }

    $self->{'data'}->{'range_low'}  = \@lows;
    $self->{'data'}->{'range_high'} = \@highs;
    
    return $self;
}




# Compute the box-count $q$-powered occupancy for a given radius.
# This is the non-optimized version.
#
# Standard usage:
# 
# $fdObj->box_count(0.2);
#
# It stores the $S_{q}(r)$ value inside the object...
sub box_put($$) {
    my $self   = shift;
    my $radius = shift;
    my $data_o = undef;
    my $dims   = undef;
    my $count  = undef;
    my $q      = undef;
    my @highs  = undef;
    my @lows   = undef; 

    my %table   = ();
    my $normalize = $self->{'params'}->{'normalize'};
    my $obj_idx = undef;


    if (!exists($self->{'data'}->{'box'})) {
        $self->{'data'}->{'box'} = +{};
    }


    my ($tfh, $tnam) = mktemp();
    my $dbh = tie %table, 'DB_File', $tnam, O_RDWR|O_TRUNC, 0700, $DB_BTREE;


    if (!exists($self->{'params'}->{'data_o'})) {
        die "Corrupted object, doesn't even have an undef data_o";
    }

    if (!exists($self->{'params'}->{'data_o'})) {
        die "Incomplete object, needs a defined data_o reference.";
    }

    $data_o = $self->{'params'}->{'data_o'};
    $q      = $self->{'params'}->{'q'};
    $dims   = $data_o->get_dims();
    $count  = $data_o->get_count();
    
    @lows   = @{$self->{'data'}->{'range_low'}};
    @highs  = @{$self->{'data'}->{'range_high'}};
    
    for ($obj_idx=0; $obj_idx < $count; $obj_idx++) {
        my @obj_data = $data_o->get_obj($obj_idx);      
        my $cell_tag = undef;
        my $dim_idx  = undef;
        my @tag_data = ();

        defined(@obj_data) || die "Failed to retrieve data obj $obj_idx";
        (scalar(@obj_data) == $dims) || die "Unexpected size for obj $obj_idx";

        if ($normalize != 0) {
            # normalize to a hypercube
            for ($dim_idx=0; $dim_idx < $dims; $dim_idx++) {
                if ($lows[$dim_idx] == $highs[$dim_idx]) {
                    $obj_data[$dim_idx] = 0;
                } else {
                    $obj_data[$dim_idx] = 
			($obj_data[$dim_idx] - $lows[$dim_idx]) /
                        ($highs[$dim_idx] - $lows[$dim_idx]);
                }
            }
        }
        
        # compute a cell tag
	for ($dim_idx=0; $dim_idx < $dims; $dim_idx++) {
	    $tag_data[$dim_idx] = floor(($obj_data[$dim_idx] 
					 - $lows[$dim_idx])
					/ $radius);
	}

        $cell_tag = join(' ', @tag_data);

	my $temp;
        if ($dbh->get($cell_tag, $temp)) {
            $table{$cell_tag} = 1;
        } else {
            $table{$cell_tag}++;       
	}
    }   

    $self->{'data'}->{'boxdata'} = \%table;
    $self->{'data'}->{'dbref'}   = $dbh;
    $self->{'data'}->{'tnam'}    = $tnam; 
}


# Call AFTER a box_put, to find out how many were in each cell.
#
# The box occupancy data is PRESERVED, in the 'old' slot.  This
# permits a certain optimization when the radius multiplier is an
# integer. 
sub box_count($$) {
    my $self   = shift;
    my $radius = shift;
    my $cell   = undef;
    my $occ    = undef;
    my $tref   = $self->{'data'}->{'boxdata'};
    my $dbh    = $self->{'data'}->{'dbref'};
    my $q      = $self->{'params'}->{'q'};
    my $total  = 0;

    defined($tref) || die 'Failed to actually count boxes...';

    if ($q != 1) {
	while (($cell, $occ) = each %{$tref}) {
	    my $buf;
	    if (!($dbh->get($cell, $buf))) {
		$total += pow($occ, $q);
	    }
	}
    } else {
	# q=1 is based on entropy.  Negative entropy, actually.
	my $count = $self->{'data'}->{'count'};

	while (($cell, $occ) = each %{$tref}) {
	    my $buf;
	    
	    if (!($dbh->get($cell, $buf))) {
		my $prob = $occ / $count; 
		$total += $prob * log($prob);
	    }
	}

	$total *= 1 / log(2);
    }

    $self->{'box'}->{$radius} = $total;

    if (exists($self->{'data'}->{'old_dbref'})) {
	delete $self->{'data'}->{'old_dbref'};
	untie %{$self->{'data'}->{'old_boxdata'}};
	untemp $self->{'data'}->{'old_tnam'};
    }

    $self->{'data'}->{'old_dbref'}   = $self->{'data'}->{'dbref'};
    $self->{'data'}->{'old_boxdata'} = $self->{'data'}->{'boxdata'};
    $self->{'data'}->{'old_tnam'}    = $self->{'data'}->{'tnam'};

    delete $self->{'data'}->{'dbref'};
    delete $self->{'data'}->{'boxdata'};
    delete $self->{'data'}->{'tnam'};

    return $self;
}




# Compute the box-counts for every applicable radius.
#
# Standard usage:
# 
# $fdObj->all_boxes();
#
# It stores the log/log data...
sub all_boxes($) {
    my $self        = shift;
    my $radius_min  = $self->{'params'}->{'r_min'};
    my $radius_max  = $self->{'params'}->{'r_max'};
    my $radius_ct   = $self->{'params'}->{'r_count'};
    my $q           = $self->{'params'}->{'q'};
    my $radius_mult = 1;
    my $radius      = $radius_min;
    my $idx         = 0;
    my @log_r       = ();
    my @log_s       = ();
    my $log2        = log(2);

    if ($radius_ct > 1) {
        $radius_mult = exp(log($radius_max / $radius_min) / ($radius_ct - 1));
    } else {
	die "Trying to compute a fractal dimension with but a single radius is rather pointless.";
    }

    # Hey, it's an integer multiplier.  Take advantage of grid
    # alignment optimizations.
    if ($radius_mult == int($radius_mult)) {
	return $self->all_boxes_int();
    }

    my @radii = ();

    # Compute all radii.
    for ($idx=0; $idx < $radius_ct; $idx++, $radius *= $radius_mult) {
        push @radii, $radius;
    }
    
    my $mid     = int($radius_ct/2);
    my $too_lo  = $self->{'data'}->{'count'};              # all in separate boxes
    my $too_hi  = pow($too_lo, $self->{'params'}->{'q'});  # all in one box
    my $thresh  = 0;
    my $comp_mult= 1;
 
    if ($too_hi < $too_lo) {
	$comp_mult = -1;
    }


    $idx = $mid;

    # Starting in the middle, decrease radius until all points are
    # in individual cells and stay that way (grid alignment issues),
    # or we hit the minimum.
    find_low:  {
	$radius = $radii[$idx];

	$self->box_put($radius);
        $self->box_count($radius);

	my $s = $self->{'box'}->{$radius};

	unshift @log_r,  log($radius);
	if ($q != 1) {
	    unshift @log_s,  log($s)/$log2;
	} else {
	    unshift @log_s,  $s;
	}

	if ((((($s - $too_lo) * $comp_mult) > 0) || (($s == $too_lo) && ($thresh < 5))) && 
	    ($idx > 0)) {
	    $thresh += ($s == $too_lo);
	    $idx--;
	    redo find_low;
	}
    }

    $idx = $mid+1;

    # Increase until they're all in a single cell, or we hit the max.
    if ($idx != $radius_ct) {
	$thresh = 0;
      find_high:  {
	  $radius = $radii[$idx];

	  $self->box_put($radius);
	  $self->box_count($radius);

	  my $s = $self->{'box'}->{$radius};
	  
	  push @log_r,  log($radius);
	  if ($q != 1) {
	      push @log_s,  log($s)/$log2;
	  } else {
	      push @log_s,  $s;
	  }

	  if (((($s - $too_hi) * $comp_mult) <  0) || (($s == $too_hi) && ($thresh < 5))
	      && ($idx < ($radius_ct - 1))) {
	      $idx++;
	      $thresh += ($s == $too_hi);
	      redo find_high;
	  } 
      }
    }

    $self->{'data'}->{'log_r'} = \@log_r;
    $self->{'data'}->{'log_s'} = \@log_s;

    if (exists($self->{'data'}->{'old_dbref'})) {
	delete $self->{'data'}->{'old_dbref'};
	untie %{$self->{'data'}->{'old_boxdata'}};
	untemp $self->{'data'}->{'old_tnam'};
    }

    return $self;
}




# If we have an integer radius multiplier, and normalize the same
# data in the same way... and we have the previous box count data,
# we can generate the next one faster -- iterate over occupied
# cells instead of objects.  This'll work because with an integer
# multiplier, we can overlay larger grids on the existing grid
# pattern.
#
# It requires, and assumes, that the radius is the previous *
# the radius multiplier -- otherwise this optimization is
# impossible.
sub box_put_int($) {
    my $self       = shift;
    my $cell       = undef;
    my $occ        = undef;
    my $old_tref   = $self->{'data'}->{'old_boxdata'};
    my $old_dbh    = $self->{'data'}->{'old_dbref'};
    my $q          = $self->{'params'}->{'q'};
    my %table      = ();

    my $radius_mult = $self->{'data'}->{'radius_mult'};

    defined($old_tref) || die 'Failed to actually count boxes...';

    my ($tfh, $tnam) = mktemp();
    my $dbh = tie %table, 'DB_File', $tnam, O_RDWR|O_TRUNC, 0700, $DB_BTREE;

    while (($cell, $occ) = each %{$old_tref}) {
	my $buf;

	if (!($old_dbh->get($cell, $buf))) {
	    # occupied cell

	    my @cell_idx = split / /, $cell;
	    
	    # compute new cell tag
	    @cell_idx = map { floor($_ / $radius_mult) } @cell_idx;

	    my $newbuf;
	    my $newcell = join ' ', @cell_idx;

	    if ($dbh->get($newcell, $newbuf)) {
		$table{$newcell}  = $buf;
	    } else {
		$table{$newcell} += $buf;
	    }
	}
    }

    $self->{'data'}->{'boxdata'} = \%table;
    $self->{'data'}->{'dbref'}   = $dbh;
    $self->{'data'}->{'tnam'}    = $tnam;  

    return $self;
}




# This version is ONLY meant for those with integral radius_mult.
sub all_boxes_int($) {
    my $self        = shift;
    my $radius_min  = $self->{'params'}->{'r_min'};
    my $radius_max  = $self->{'params'}->{'r_max'};
    my $radius_ct   = $self->{'params'}->{'r_count'};
    my $q           = $self->{'params'}->{'q'};
    my $radius_mult = 1;
    my @log_r       = ();
    my @log_s       = ();
    my $log2        = log(2);

    if ($radius_ct > 1) {
        $radius_mult = exp(log($radius_max / $radius_min) / ($radius_ct - 1));
    }

    my $radius_mid  = $radius_min * ($radius_mult ** (int($radius_ct / 2))); 
    my $radius      = $radius_mid;

    ($radius_mult == int($radius_mult)) || 
	die 'all_boxes_int called with a float value';

    $self->box_put($radius);
    $self->box_count($radius);

    my $s  = $self->{'box'}->{$radius};
    my $ct = $self->{'data'}->{'count'};

    push @log_r,  log($radius)/$log2;

    if ($q != 1) {
	push @log_s,  log($s)/$log2;
    } else {
	push @log_s,  $s;
    }

    my $max_bc = pow(($self->{'data'}->{'count'}), $self->{'params'}->{'q'});
    my $min_bc = $self->{'data'}->{'count'};
    my $comp_mult = ($max_bc > $min_bc) ? 1 : -1;
    
    $self->{'data'}->{'radius_mult'} = $radius_mult;

    # Increase radius until we hit the max, or it'd be pointless
    # (all objects in same cell).
    for ($radius=$radius * $radius_mult; $radius <= $radius_max; 
	 $radius *= $radius_mult) {

	# box_put_int(), not box_put(); _int is faster going
	# forwards.
	$self->box_put_int();
	$self->box_count($radius);

	$s = $self->{'box'}->{$radius};
	push @log_r,  log($radius)/$log2;

	if ($q != 1) {
	    push @log_s,  log($s)/$log2;
	} else {
	    push @log_s,  $s;
	}

	if ($s == $max_bc) {
	    last;
	}
    }

    # Likewise, reduce.
    if ((($self->{'box'}->{$radius_mid} - $min_bc) * $comp_mult) > 0) {
	for ($radius=$radius_mid / $radius_mult; $radius >= $radius_min;
	     $radius /= $radius_mult) {
	    $self->box_put($radius);
	    $self->box_count($radius);


	    $s = $self->{'box'}->{$radius};
	    unshift @log_r,  log($radius)/$log2;
	    if ($q != 1) {
		unshift @log_s,  log($s)/$log2;
	    } else {
		unshift @log_s, $s;
	    }
	    if ($s == $min_bc) {
		last;
	    }
	}
    }

    $self->{'data'}->{'log_r'} = \@log_r;
    $self->{'data'}->{'log_s'} = \@log_s;

    if (exists($self->{'data'}->{'old_dbref'})) {
	delete $self->{'data'}->{'old_dbref'};
	untie %{$self->{'data'}->{'old_boxdata'}};
	untemp $self->{'data'}->{'old_tnam'};
    }
    
    return $self;
}



# If the packages is invoked with  
#
#   method => 'slow'
#
# we compute all pairwise distances and sort them to generate
# log/log data.  This is both slow and memory-intensive, as we
# no longer rely on temporary files, and will end up doing an
# in-memory sort of all pairwise distances. 
#
# The parameter $q is ignored; this is aimed soley at equivalence
# to the fractal dimension.
#
# Straight Euclidean distance is used.
#
# Just to emphasize:
#
#     **** *    **** *   *  
#     *    *    *  * *   *
#     **** *    *  * * * * 
#        * *    *  * * * *
#     **** **** ****  * *
#
# DO NOT USE on anything except small data sets, unless you have
# a lot of time and memory.

sub all_pairs($) {
    my $self      = shift;
    my $ct        = $self->{'data'}->{'count'};
    my $dims      = $self->{'data'}->{'dims'};
    my $data_obj  = $self->{'params'}->{'data_o'};
    my @distances = ();
    my $i         = 0;
    my $zero    = 0;
    my $log2    = log(2);

    for ($i=0; $i < ($ct-1); $i++) {
	my @i_obj = $data_obj->get_obj($i);
	my $j     = $i+1;

	for ($j=$i+1; $j < $ct; $j++) {
	    my $dist  = 0;
	    my @j_obj = $data_obj->get_obj($j);
	    my $k     = 0;

	    for ($k=0; $k < $dims; $k++) {
		$dist += ($i_obj[$k] - $j_obj[$k]) ** 2;
	    }	    

	    $dist = sqrt($dist);

	    if ($dist == 0) {
		$zero++;
	    } else {
		push @distances, $dist;
	    }
	}
    }

    @distances = sort { $a <=> $b } @distances;

    my $cumulative = $zero;
    my $dist       = 0;

    my $log_r_ref = ($self->{'data'}->{'log_r'} = +[]);
    my $log_s_ref = ($self->{'data'}->{'log_s'} = +[]);

    distloop:  {
	if ((!defined(@distances)) || (scalar(@distances) < 1)) {
	    last distloop;
	}

	my $this_dist = shift @distances;
	$cumulative++;

      matchloop:  {
	  if ((!defined(@distances)) || (scalar(@distances) < 1)) {
	      last matchloop;
	  }

	  if ($distances[0] == $this_dist) {
	      $cumulative++;
	      shift @distances;
	      redo matchloop;
	  }

	  push @$log_r_ref, log($this_dist)/$log2;
	  push @$log_s_ref, log($cumulative)/$log2;
      }

	redo distloop;	
    }

    return $self;
}




# If the method is 'slow' (namely, pairwise evaluation), we expect
# a certain density bias towards the upper end.  This function 
# attempts to alleviate that by dividing the data into intervals
# according to log(pairwise distance), and picking at most one
# point from each interval.

sub filter_points($) {
    my $self  = shift;
    my @log_r = @{$self->{'data'}->{'log_r'}};
    my @log_s = @{$self->{'data'}->{'log_s'}};
    my $intervals = $self->{'params'}->{'intervals'};

    defined($intervals) || die "Shouldn't be here without defining param.";

    (defined(@log_r) && (scalar(@log_r) > 0)) || die "Sanity check failed!";
    (defined(@log_s) && (scalar(@log_s) == (scalar(@log_r)))) || 
     die "Sanity check failed!";

    my $min_dist  = $log_r[0];
    my $max_dist  = $log_r[-1];
    my @new_log_r = ($log_r[0]);
    my @new_log_s = ($log_s[0]);
    my $interval  = 1;
    my $idx       = 1;
    my $max       = scalar(@log_r);


    if ($intervals < 2) {
	# Why are you doing this???  WHY??
	
	$self->{'data'}->{'log_r'} = \@new_log_r;
	$self->{'data'}->{'log_s'} = \@new_log_s;
	return $self;
    }

    my $int_step  = ($max_dist - $min_dist) / $intervals;
    my $int_min   = $int_step + $min_dist;

    for ($interval=1, $idx=1; ($interval < $intervals) &&
	 ($idx < $max); $idx++) {

	my $log_r = $log_r[$idx];
	my $log_s = $log_s[$idx];

	if ($log_r >= $int_min) {
	    push @new_log_r, $log_r;
	    push @new_log_s, $log_s;

	  find_next_interval:  {
	      $int_min += $int_step;
	      $interval++;
	      if ($log_r >= $int_min) {
		  redo find_next_interval;
	      }
	  }
	} 
    }

    $self->{'data'}->{'log_r'} = \@new_log_r;
    $self->{'data'}->{'log_s'} = \@new_log_s;
    return $self;
}


# Trim the flat portions at the beginning and the end of the
# log_r and log_s arrays.
#
# Standard usage:
# 
# $fdObj->trim_ranges();
#
# It stores the log/log data...
sub trim_ranges($) {
  my $self  = shift;
  my @log_r = @{$self->{'data'}->{'log_r'}};
  my @log_s = @{$self->{'data'}->{'log_s'}};
  my $diff  = 0.01;
  my $min   = $log_s[0];
  my $max   = $log_s[0];
  my $i     = 0;
  my $j     = scalar @log_s;

  if ($j < 5) {
      return;
  }

  for ($i=1; $i < $j; $i++) {
      $min = ($min > $log_s[$i]) ? $log_s[$i] : $min;
      $max = ($max < $log_s[$i]) ? $log_s[$i] : $max;
  }
  
  $min = $min+$diff;
  $max = $max-$diff;

  my @new_log_r = ();
  my @new_log_s = ();
  my $line_ct   = grep { ($_ >= $min) && ($_ <= $max) } @log_s;

  if ($line_ct > 5) {
      for ($i=0; $i < $j; $i++) {
	  if ($log_s[$i] < $min) {
	      next;
	  }
	  
	  if ($log_s[$i] > $max) {
	      next;
	  }

	  push @new_log_r, $log_r[$i];
	  push @new_log_s, $log_s[$i];
      }
  } else {
      for ($i=0; $i < $j; $i++) {
	  if (($i < ($j-1)) && ($log_s[$i] == $log_s[$i+1])) {
	      next;
	  }

	  push @new_log_r, $log_r[$i];
	  push @new_log_s, $log_s[$i];
      }
  }

  $self->{'data'}->{'log_r'} = \@new_log_r;
  $self->{'data'}->{'log_s'} = \@new_log_s;

  return $self;
}


# report log/log points (the final ones fed next into line-fitting)
sub print_loglog($) {
    my $self  = shift;
    my @log_r = @{$self->{'data'}->{'log_r'}};
    my @log_s = @{$self->{'data'}->{'log_s'}};
    my $r_ct  = scalar @log_r;
    my $s_ct  = scalar @log_s;
    my $idx   = 0;

    ((defined($r_ct)) && (defined($s_ct)) && ($r_ct == $s_ct)) ||
	die "Strange log-log data!";

    print ";; ===== begin-log-log =====\n";
    
    for ($idx=0; $idx < $r_ct; $idx++) {
	my $log_r = $log_r[$idx];
	my $log_s = $log_s[$idx];

	print "; $log_r $log_s\n";
    }
    print ";; ===== end-log-log ===== \n";
}


# Taken from Christos's original fractal dimension package,
# with minor modifications for usage here.
# 
# Returns (slope, correlation coefficient)
sub lsfit($$$) {
  my $self  = shift;
  my $xref  = shift;
  my $yref  = shift;
  my $x1    = undef;    # $\sum {x}$
  my $x2    = undef;    # $\sum {x}^{2}$
  my $y1    = undef;    # $\sum {y}$
  my $y2    = undef;    # $\sum {y}^{2}$
  my $xy    = undef;    # $\sum ({x} \times{y})$
  my $n     = undef;    # count
  my $i     = undef;    # index
  my $xlen  = scalar @{$xref};
  my $ylen  = scalar @{$yref};
  my $xval  = undef; 
  my $yval  = undef;
  
  my $a     = undef;
  my $b     = undef; 
  my $r     = undef;   # y = ax+b
  
  
  if ($xlen != $ylen) {
    die "lsfit:  log_r is of length '$xlen', log_s '$ylen'";
  }
  
  $n = $xlen;
  
  $x1 = $x2 = $y1 = $y2 = $xy = 0;
  
  for($i=0; $i<$n; $i++) {
    $xval = $xref->[$i];
    $yval = $yref->[$i];
    $x1 += $xval;
    $y1 += $yval;
    $x2 += $xval * $xval;
    $y2 += $yval * $yval;
    $xy += $xval * $yval;
  }
  
  if ($n > 1) {
    $a = ($xy - $x1 * $y1 / $n) / ($x2 - $x1 * $x1 / $n);
    $b = ($y1 - $a * $x1) / $n;
    
    # correlation coefficient r
    my $xvar = sprintf '%.5f', $x2 - $x1 * $x1 / $n;
    my $yvar = sprintf '%.5f', $y2 - $y1 * $y1 / $n;
    if($xvar == 0) { 
      $r = 1;
    } elsif( $yvar == 0 ){ 
      $r = 1;
    } else {
      $r = ($xy - $x1 * $y1 / $n) / sqrt($xvar) / sqrt($yvar);
    }
    
    # old version:
    # $r = ($xy - $x1 * $y1 / $n) / sqrt($x2 - $x1 * $x1 / $n) 
    # / sqrt($y2 - $y1 * $y1 / $n);
    # print 'slope= ', $a, '    y_intcpt= ', $b, '  corr= ', $r;
  } else {
    # print 'slope= ', 9999, '  y_intcpt= ', 0, '   corr= ', 0;
    $a = 9999; $b=0; $r=0;
  }
  
  return( $a, $b, $r);
}                               # end lsfit



sub lts($) {
    my $self  = shift;
    my $xaref = $self->{'data'}->{'log_r'};   # x values
    my $yaref = $self->{'data'}->{'log_s'};   # y values
    my ($n, $i);
    my $xlen;
    my $ylen;
    my (@smallx, @smally);
    my ($a, $b, $r);   # intercept, slope and corr. coeff.
    my ($start, $end); # offsets of matching portion
    my $q;	# stretch to match = end-start +1
    
    my $diff;
    my $currentDiff ;
    my $currentStart ;
    my ( $currenta , $currentb , $currentr);
    
    
    $xlen = scalar( @$xaref);
    $ylen = scalar( @$yaref);
    
    $n = $xlen;
    
    $q = int ( $n/2) +1; #according to Venable + Ripley
    #$q = int($n/2);

    for( $start =0; $start < $n - $q + 1 ; $start ++ ){
	$end = $start + $q -1;
	
	# also works
	# @smallx = @{ $xaref} [ $start .. $end ];
	# @smally = @{ $yaref} [ $start .. $end ];
	
	my $i;
	for($i=0; $i<$q; $i++){
	    $smallx[$i] = $$xaref[$i+$start];
	    $smally[$i] = $$yaref[$i+$start];
	}
	
	($a, $b, $r) = $self->lsfit(\@smallx , \@smally );
	# surprisingly, \@smallx etc DOES NOT WORK!
	# $a=$b=$r= 9999999; #TST
	
	
	$diff = sumsqdiff( \@smallx, \@smally, $a, $b);
	if( !defined($currentDiff) || ($diff < $currentDiff ) ) {
	    # new champion, with smaller diffs
	    $currentDiff = $diff;
	    $currentStart = $start;
	    $currenta = $a;
	    $currentb = $b;
	    $currentr = $r;
	}
	
    }# end for
	
	return( $currenta, 
	       $currentb,
	       $currentr);

}# end lts




# given the matrix @x and @y
# and the fitting line y=ax+b
# it computes the sum of squared differences
#    sum ( y_i - a*x_i - b ) ^ 2
sub sumsqdiff($$$$){
    my $xaref = shift;
    my $yaref = shift;
    my $a     = shift;
    my $b     = shift;
    my $xlen;
    my $ylen;
    my $res;

    # parse the parms
    $xlen = scalar( @$xaref);
    $ylen = scalar( @$yaref);

    die unless ($xlen == $ylen);

    my $n = $xlen;
    my $i;
    $res =0;
    for($i=0; $i<$n; $i++){
	my $xval = $$xaref[$i];
	my $yval = $$yaref[$i];
        $res += ($yval - $a* $xval - $b) ** 2;
    }

    return($res);
}




