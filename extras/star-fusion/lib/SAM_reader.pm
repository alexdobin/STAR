package SAM_reader;

use strict;
use warnings;
use Carp;

use SAM_entry;

sub new {
	my $packagename = shift;
	my $filename = shift; 

	unless ($filename) {
		confess "Error, need SAM filename as parameter";
	}
	
	my $self = { filename => $filename,
				 _next => undef,
				 _fh =>  undef,
	};

	bless ($self, $packagename);

	$self->_init();

	return($self);
}


####
sub _init {
	my ($self) = @_;
    
    if ($self->{filename} =~ /\.bam$/) {
        open ($self->{_fh}, "samtools view $self->{filename} |") or confess "Error, cannot open file " . $self->{filename};
    }
    elsif ($self->{filename} =~ /\.gz$/) {
        open ($self->{_fh}, "gunzip -c $self->{filename} | ") or confess "Error, cannot open file " . $self->{filename};
    }
    else {
        open ($self->{_fh}, $self->{filename}) or confess "Error, cannot open file " . $self->{filename};
    }
    
	$self->_advance();

	return;
}

####
sub _advance {
	my ($self) = @_;

	my $fh = $self->{_fh};

	my $next_line = <$fh>;
	while (defined ($next_line) && ($next_line =~ /^\@/ || $next_line !~ /\w/)) {  ## skip over sam headers
		$next_line = <$fh>;
	}
	
	if ($next_line) {
		$self->{_next} = new SAM_entry($next_line);
	}
	else {
		$self->{_next} = undef;
	}
	
	return;
}

####
sub has_next {
	my $self = shift;
	
	if (defined $self->{_next}) {
		return(1);
	}
	else {
		return(0);
	}
}


####
sub get_next {
	my $self = shift;
	
	my $next_entry = $self->{_next};

	$self->_advance();
	
	if (defined $next_entry) {
		return($next_entry);
	}
	else {
		return(undef);
	}
}

####
sub preview_next {
	my $self = shift;
	return($self->{_next});
}


1;

	
	
