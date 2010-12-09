#!/usr/bin/perl -w

my $idsfile = $ARGV[0];
my $readsfile = $ARGV[1];

open ID, $idsfile;
open READ, $readsfile;


while(my $id = <ID>) {
    chomp $id;
    print ">$id\n";
    my $read = <READ>;
    chomp $read;
    $read =~ s/^\s+/'-' x (length($&))/e;
    print $read,"\n";  
}
