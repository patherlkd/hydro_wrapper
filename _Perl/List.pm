#!/usr/bin/perl -w

package _Perl::List;

use v5.10;
use strict;
use Exporter qw{ import };

our @EXPORT=qw( contains );

sub contains{

    my ($ref,$element) = @_;

    my @array = @{$ref};
    my $bool= 0;

    for my $x (@array){
        if($x eq $element){
            $bool= 1;
        }
    }
    return $bool;
}

1;
