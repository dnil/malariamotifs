#!/usr/bin/perl -w 


my $n_co = 61;
my $n_se = 41;
my $n_seq = 91;
my $dominance_tab = "Rosetting.dominance060313.1dom.tab";

my $n_ass = $n_co + $n_se;

my @non_random_assignments;

push @non_random_assignments, split (/ */,'S'x$n_se);
push @non_random_assignments, split (/ */,'C'x$n_co);

my @new_order; 

for(my $i = 0; $i < $n_ass ; $i++) {
    $random_element = rand($n_ass - $i) ;

    push @new_order, splice @non_random_assignments, $random_element, 1;
}

# all slots guaranteed to be filled once. later randomly distribute remainder.

my @seq_assignment;

for (my $i = 0; $i < $n_seq; $i++) {
    if($new_order[$i] eq 'S') {
	$seq_assignment[$i] = "se00";
    } elsif($new_order[$i] eq 'C') {
	$seq_assignment[$i] = "co00";
    }
}

for(my $i= $n_seq; $i < $n_ass; $i++) {
    $add_to_sequence_nr = rand($n_seq);

    if($new_order[$i] eq 'S') {
	$seq_assignment[$add_to_sequence_nr] .= "_se00";
    } elsif($new_order[$i] eq 'C') {
	$seq_assignment[$add_to_sequence_nr] .= "_co00";
    }    
}

open TAB, $dominance_tab;

while(<TAB>) {
    my ($seq_name) = m/^(\S+)\s+\S+/;
    print $seq_name, "\t", shift @seq_assignment, "\n";
}
