: # Use perl...
eval 'exec perl5 -w -S $0 "$@"' 
if 0;

#!/usr/bin/perl5 -w
#
############################################################
# $Id: DiskFracTest.pl,v 1.3 2000/10/24 18:00:46 lw2j Exp $
# 
# $Log:	DiskFracTest.pl,v $
# Revision 1.3  2000/10/24  18:00:46  lw2j
# Standardized path.
# 
# Revision 1.2  2000/10/03  14:08:01  lw2j
# Added magic perl invoker.
# 
# Revision 1.1  2000/07/18  13:27:11  lw2j
# Initial revision
#   
############################################################
#
# A simple demo to invoke the FracDim package on a single file.
# It's not terribly nice coding.
#
use strict;

##### <change_me>
#
# Change this to where you placed the modules.
#
use lib '/usr0/lw2j/private/Research/FracDim';
##### </change_me>

use Symbol;

# DiskWrapper can actually be substituted for, actually, if 
# you've got a module with a similar interface.

require DiskWrapper;
require DiskFracDim;

my $tri_fh  = gensym;
my $obj     = undef;
my $fd_obj  = undef;
my %params  = ();
my $fname   = undef;


# parse arguments
{
    my $arg;

    while (defined($arg = $ARGV[0]) && ($arg =~ /^-/)) {
	$arg = shift @ARGV;

	if ($arg eq '-v') {
	    $params{'verbose'} = 1;
	} elsif ($arg eq '-s') {
	    $params{'method'}  = 'slow';
	} elsif (($arg eq '-q')       ||
		 ($arg eq '-r_min')   ||
		 ($arg eq '-r_max')   ||
		 ($arg eq '-r_count') ||
		 ($arg eq '-intervals')) {
	    
	    $arg =~ s/^-//;
	    my $val = shift @ARGV;
	    
	    defined($val) || die "Failed to supply value for -$arg";

	    if ($val ne 'undef') {
		$params{$arg} = $val; 
	    } else {
		$params{$arg} = undef;
	    }
	} elsif ($arg eq '--') {
	    last;
	} else {
	    print_usage_and_exit();
	}
    }
    $fname   = shift @ARGV;
}

if (!defined($fname)) {
    print_usage_and_exit();
}

open($tri_fh, $fname) || die "Failed to open '$fname':  $!";
$obj = DiskWrapper::new DiskWrapper($tri_fh);
$params{'data_o'} = $obj;
$fd_obj = new DiskFracDim(%params);
my ($fd, $y_int, $corr) = $fd_obj->fracdim();

print "fd = $fd, y_int = $y_int, corr = $corr\n";
close($tri_fh);

exit 0;


#####################
sub print_usage_and_exit() {
    use File::Basename;

    my $progname = basename($0);

    print "$progname -- Interact with the DiskFracDim.pm module\n";
    print "Usage:\n";
    print "\n";
    print "$progname [-v] [-s] [-q nn] [-r_min nn] [-r_max nn]\n";
    print "          [-r_count] [-intervals]  filename\n";
    print "\n";
    print "Exact parameter semantics are defined by the module.\n";
    print "\n";
    print " -v          Specifies 'verbose'.\n";
    print " -s          Use the slow all-pairs method for C2 -- more\n";
    print "             precise, but perhaps unbearably slow and vile\n";
    print "             in terms of memory consumption.\n";
    print " -q          Box-counting exponent.\n";
    print " -r_min      Minimum radius the module will consider using.\n";
    print " -r_max      Maximum radius it'll try.\n";
    print " -r_count    Maximum number of radii it'll try.\n";
    print " -intervals  The number of intervals for the all-pairs method.\n";
    print "\n";
    print "The options after -s are all specific to box-counting.\n";
    print "\n";
    print "The filename must correspond to a vector data file, with\n";
    print "one equal-dimensioned vector per line, separated by\n";
    print "whitespace or commas.  Comments are to be delineated with\n";
    print "semicolons.\n";
    exit 0;
}
