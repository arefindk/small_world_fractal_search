: # Find perl...
eval 'exec perl5 -w -S $0 "$@"' 
if 0;

#!/usr/bin/perl5 -d
#
############################################################
# $Id: mktemp.pm,v 1.2 2000/10/03 14:08:01 lw2j Exp $
# 
# $Log:	mktemp.pm,v $
# Revision 1.2  2000/10/03  14:08:01  lw2j
# Added magic perl invoker.
# 
# Revision 1.1  2000/07/18  13:27:11  lw2j
# Initial revision
#   
############################################################
#
# Purpose:  
#  To supply a mktemp()-like function.
#
##################################################################
#
# Package Header
#
##################################################################

package mktemp;
require Exporter;

##########################################################################
my %mktemp_tmpfiles;


BEGIN {
	use Fcntl;
	use POSIX qw(tmpnam);	
	use Symbol;
	
	use Exporter ();
	use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

	$VERSION     = 1.0;
	@ISA         = qw(Exporter);
	@EXPORT      = qw(mktemp untemp untempAll forget);
	%EXPORT_TAGS = ();
	@EXPORT_OK   = ();


	%mktemp_tmpfiles = ();

	return TRUE;
 }

##########################################################################


# Erase a file created with mktemp.
sub untemp($) {
	my($name) = @_;
	
	if (!exists($mktemp_tmpfiles{$name})) {
		die "untemp invoked with non-mktemp file '$name'";
	}
	
	unlink($name) || die "Failed to unlink '$name': $!";
	delete $mktemp_tmpfiles{$name};
}


# Pretend we never made the file -- perhaps the programmer wants
# to keep it around for a bit longer, and doesn't want it to be
# auto-removed.
sub forget($) {
	my($name) = @_;

	delete $mktemp_tmpfiles{$name};
}


# Calls untemp on all temp files we didn't forget().
sub untempAll() {
	my $name;
	
	while (defined($name = (each %mktemp_tmpfiles))) { 
		&untemp($name);
	}
}



# Make a temporary file -- tmpnam() generates in /tmp, at least
# on Unices.
sub mktemp() {
	my($fh);
	my($name);

	# Minor precaution.  Remove if you prefer different
	# sig handlers.
	if (!defined($SIG{'HUP'})) {
		$SIG{'HUP'} = 'untempAll';
	}
	
	if (!defined($SIG{'KILL'})) {
		$SIG{'KILL'} = 'untempAll';
	}
	
	$fh = &gensym();

	# Be stubborn.
	do { 
	    $name = &tmpnam();
	} until (sysopen($fh, $name, O_RDWR|O_CREAT|O_EXCL));
	
	$mktemp_tmpfiles{$name} = 1;
	
	return($fh, $name);
}




END 
	{ 
		&untempAll();
	}

##########################################################################

return TRUE;
