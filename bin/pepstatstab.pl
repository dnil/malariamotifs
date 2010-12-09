#!/usr/bin/perl -w

while(<STDIN>) {
    if (m/Residues = (\d+)/) { print "Residues\t$1\n"; }
    elsif (m/=\s+(\w+)\s+\d+\s+([\d+\.]+)\s+\d+/) { print $1,".moleperc\t", $2,"\n"; }
    elsif (m/^(\w+)\s+.+\d+\s+([\d+\.]+)/) { print "$1.moleperc\t$2\n"; }

    if(m/^C = Cys\s+(\d+)/) { print "Cyscount\t$1\n"; }
}
