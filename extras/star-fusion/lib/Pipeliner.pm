package Pipeliner;

use strict;
use warnings;
use Carp;

################################
## Verbose levels:
## 1: see CMD string
## 2: see stderr during process
my $VERBOSE = 0;
################################


####
sub new {
    my $packagename = shift;
    my %params = @_;
    
    if ($params{-verbose}) {
        $VERBOSE = $params{-verbose};
    }
    
    my $self = { 
        cmd_objs => [],
    };
    
    bless ($self, $packagename);

    return($self);
}


sub add_commands {
    my $self = shift;
    my @cmds = @_;
    
    foreach my $cmd (@cmds) {
        unless (ref($cmd) =~ /Command/) {
            confess "Error, need Command object as param";
        }
        push (@{$self->{cmd_objs}}, $cmd);
    }
    
    return $self;

}

sub run {
    my $self = shift;

    foreach my $cmd_obj ($self->_get_commands()) {
        
        my $cmdstr = $cmd_obj->get_cmdstr();
        my $checkpoint_file = $cmd_obj->get_checkpoint_file();

        if (-e $checkpoint_file) {
            print STDERR "-- Skipping CMD: $cmdstr, checkpoint exists.\n" if $VERBOSE;
        }
        else {
            print STDERR "* Running CMD: $cmdstr\n" if $VERBOSE;
            
            my $tmp_stderr = "tmp.$$.stderr";
            if (-e $tmp_stderr) {
                unlink($tmp_stderr);
            }
            unless ($VERBOSE == 2) {
                $cmdstr .= " 2>$tmp_stderr";
            }
            
            my $ret = system($cmdstr);
            if ($ret) {
                
                if (-e $tmp_stderr) {
                    system("cat $tmp_stderr");
                    unlink($tmp_stderr);
                }
                                
                confess "Error, cmd: $cmdstr died with ret $ret";
            }
            else {
                `touch $checkpoint_file`;
                if ($?) {
                    
                    confess "Error creating checkpoint file: $checkpoint_file";
                }
            }

            if (-e $tmp_stderr) {
                unlink($tmp_stderr);
            }
        }
    }

    return;
}

sub _get_commands {
    my $self = shift;

    return(@{$self->{cmd_objs}});
}

package Command;
use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    
    my ($cmdstr, $checkpoint_file) = @_;

    unless ($cmdstr && $checkpoint_file) {
        confess "Error, need cmdstr and checkpoint filename as params";
    }

    my $self = { cmdstr => $cmdstr,
                 checkpoint_file => $checkpoint_file,
    };

    bless ($self, $packagename);

    return($self);
}
    
####
sub get_cmdstr {
    my $self = shift;
    return($self->{cmdstr});
}

####
sub get_checkpoint_file {
    my $self = shift;
    return($self->{checkpoint_file});
}



1; #EOM
