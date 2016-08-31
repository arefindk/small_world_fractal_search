: #Find perl...
eval 'exec perl5 -w -S $0 "$@"' 
if 0;

#!/usr/bin/perl5 -w
#
############################################################
# $Id: DiskWrapper.pm,v 1.2 2000/10/03 14:08:01 lw2j Exp $
# 
# $Log:	DiskWrapper.pm,v $
# Revision 1.2  2000/10/03  14:08:01  lw2j
# Added magic perl invoker.
# 
# Revision 1.1  2000/07/18  13:27:11  lw2j
# Initial revision
#   
############################################################
#
# Purpose:  
#  To provide a very basic wrapper for small data files.  
#
#  The data is expected to be in the form of vectors, with one
#  vector per line and all vectors being of equal dimensionality.
#  A valid filehandle is required, as the data is *not* stored
#  in memory (well, a few vectors at a time are) in order to
#  minimize memory consumption.  This filehandle must be 
#  seek-able (backwards).  
# 
#  Data may be separated by either whitespace or commas.
#  Beginning/trailing whitespace will be ignored.
#  Semicolons and octothorpes may be used to denote comments.
#
#     get_count($)       # return how many items there are
#     get_dims($)        # return embedding dimensionality
#     get_obj($$)        # return the object by 0-based index, as an
#                        # array of scalars
#
#     set_mask(@)        # The 'mask' is a list of dimensions to
#                        # be used -- if specified.  If not, then
#                        # all are used.  If the mask exists, the
#                        # attributes will basically not exist for
#                        # any purpose until the mask is reset.
#     get_mask(@)        # Returns the current mask.
#     del_mask(@)        # Deletes the current mask.
# 
#  The primary purpose of this wrapper is to provide a simple 
#  interface for, say, the Fractal Dimension package, so it should
#  be trivial to imitate the same interface for your particular
#  data needs.  If you want, you can use the masking functionality
#  to ignore given attributes -- it's a crude implementation in
#  that one could enhance it to allow for re-ordering and 
#  duplication, as well.
#
##########################################################################
#
# Package Header
#
##########################################################################

package DiskWrapper;
require Exporter;

use strict;
use POSIX;

##########################################################################

BEGIN {
    use Exporter ();
    use Symbol;
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

    $VERSION     = 1.0;
    @ISA         = qw(Exporter);
    @EXPORT =  qw(
		  new
		  load_file

		  get_count
		  get_dims
		  get_obj
		  set_mask

		  get_mask
		  del_mask
		  );
    %EXPORT_TAGS = ();
    @EXPORT_OK   = ();

    return 1;
}



# new:  creates blessed DiskWrapper objects with default params
# Supply a file reference, please.
sub new($@) {
    my $class    = shift;
    my $file_ref = shift;
    my $self     = +{};

# For instance.  The cacheing code is actually pretty irrelevant
# for most usage -- especially sequential access.  It's also the
# most trivial 'algorithm' possible.

    my $csize    = 1024;
    bless $self, $class;

    $self->{"data"}   = +{};
    $self->{"params"} = +{};
    $self->{"lineno"} = 0;
    $self->{"cachesize"} = $csize;
    $self->{"cache"}     = +[(undef) x $csize];
    
    if (defined($file_ref)) {
        $self->load_file($file_ref);
    }
    
    return $self;   
}



# Load the file.  There's some minor munging to handle both
# whitespace and commas.
sub load_file($$) {
    my $self     = shift;
    my $file_ref = shift;
    my $dims     = undef;
    my $count    = 0;
    my $line     = undef;
    my @obj      = undef;

    $self->{"data"} = +{};
    $self->{"fh"}   = $file_ref;

    while (defined($line = <$file_ref>)) {

	# This bit o' nastiness should handle most numerical ASCII 
	# flat-files.
	
        $line =~ s/,/ /go;
        $line =~ s/^\s+//o;
        $line =~ s/\s\s+/ /go;
	$line =~ s/\;.*$//go;
	$line =~ s/\#.*$//go;
        chomp($line);

        (defined($line) && ($line ne "")) || next;

        @obj = split /\s+/, $line;
#        $self->{"data"}->{$count} = +[ @obj ];
        $count++;

        if (!defined($dims)) {
            $dims = scalar @obj;
        } else {
            ($dims == (scalar @obj)) || die "Irregularly sized objects: '$line'";
        }
    }

# If it reached this point, we know that every vector is of equal length.
    $self->{"params"}->{"count"} = $count;
    $self->{"params"}->{"dims"}  = $dims;
    
    # new:  reset file pointer
    (seek $file_ref, 0, 0) || die "Failed to seek to top:  $!";
    return $self;
}




# Return the stored count... assumes the data has already been
# loaded.
sub get_count($) {
    my $self = shift;

    defined($self->{"params"}->{"count"}) || die "Can't get_count with 0 data";

    return $self->{"params"}->{"count"};
}



# Return the embedding dimensionality... assumes the data has already been
# loaded.
sub get_dims($) {
    my $self   = shift;
    my $maskct = $self->{"maskct"};

    if (defined($maskct)) {
	return $maskct;
    }

    defined($self->{"params"}->{"dims"}) || die "Can't get_dims with 0 data";

    return $self->{"params"}->{"dims"};
}



# Return a specific object, by index.
sub get_obj($$) {
    my $self  = shift;
    my $index = shift;
    my $mask_ct = $self->{"maskct"};
    
#    defined($self->{"data"}->{$index}) || die "Unexpected index '$index'";

    my $cacheidx = $index % ($self->{"cachesize"});
    my @orig     = ();
    my $fh   = $self->{"fh"};

    if (defined($self->{"cache"}->[$cacheidx])) {
	my $tag = $self->{"cache"}->[$cacheidx]->[0];
	my @data = @{$self->{"cache"}->[$cacheidx]->[1]};

	if ($tag == $index) {
	    @orig = @data;
	}
    }

    if (scalar(@orig) == 0) {
	# not in cache;
	if ($index < $self->{"lineno"}) {
	    # seek
	    (seek $self->{"fh"}, 0, 0) || die "Failed to seek to top:  $!";

	    $self->{"lineno"}=0;
	}

	cycle:  {
	    my $line = undef;

	    if (!defined($line = <$fh>)) {
		die "Failed to read to $index\n";
	    }

	    $line =~ s/\;.*$//go;
	    $line =~ s/\#.*$//go;

	    chomp($line);
	    
	    if (!($line =~ /\S/o)) {
		redo;
	    }

	    if ($self->{"lineno"} == $index) {
		$line =~ s/,/ /go;
		$line =~ s/\s+/ /go;
		@orig = split /\s/, $line;
		$self->{"lineno"}++;
		last cycle;
	    }
	    $self->{"lineno"}++;
	    redo cycle;
	}
    }

    $self->{"fh"} = $fh;
    $self->{"cache"}->[$cacheidx] = +[ $index, +[@orig] ];

    if (defined($mask_ct)) {
	my @values = ();
	my @mask   = @{$self->{"mask"}};
	my $attr;
	my $val;

	foreach $attr (@mask) {
	    push @values, $orig[$attr];
	}
	return @values;
    } else {
	return @orig;
    }
}


# Set the attribute mask.
# An undef value clears it.
sub set_mask($@) {
    my $self = shift;
    my @mask = @_;

    if (defined(@mask) && (scalar(@mask) > 0)) {
	my %table = ();
	my $attr;

	# guarantee uniqueness

	foreach $attr (@mask) {
	    $table{$attr} = 1;
	}

	@mask = sort (keys %table);
	
	$self->{"mask"}   = +[ sort(@mask) ];
	$self->{"maskct"} = scalar @mask;
    } else {
	delete $self->{"mask"};
	delete $self->{"maskct"};
    }
}              


# Return the current mask.
sub get_mask($) {
    my $self = shift;
    
    if (exists($self->{"mask"})) {
	return @{$self->{"mask"}};
    } else {
	return undef;
    }
}


# Clear the mask.
sub del_mask($) {
    my $self = shift;

    delete $self->{"mask"};
    delete $self->{"maskct"};
}
