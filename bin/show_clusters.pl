#!/usr/bin/perl -w

use Bio::TreeIO;

my $treeio = new Bio::TreeIO(-file => $ARGV[0], -format => 'newick');
my $tree = $treeio->next_tree;

my $branch_treshold = 0.015;
defined($ARGV[1]) && ($branch_treshold = $ARGV[1]);

my $rootnode = $tree->get_root_node();

my %seen_ids;

foreach my $node ($rootnode->get_all_Descendents) {

    if( $node->is_Leaf ) {
	if( $node->branch_length < $branch_treshold ) {

	    if( defined($seen_ids{$node->id}) && $seen_ids{$node->id} == 1 ) {
#	    print "Found short branch ", $node->id," bl ",$node->branch_length,"\n";
#		print "..but have already printed its companions.\n";
		next;
	    }

	    print "Found short branch ", $node->id," bl ",$node->branch_length,"\n";

	    my $in_bush = 1;
	    my ($last_parent, $next_parent);
	    $next_parent = $node;
	    while($in_bush) {
		$last_parent = $next_parent;
		$next_parent = $node->ancestor;
#		print "Get next parent ($next_parent, $last_parent, $node, $rootnode)..\n";
		if( ($next_parent == $last_parent)  || $next_parent->branch_length > $branch_treshold) {
		    $in_bush = 0;
		    print "New bush: last parent bl ",$next_parent->branch_length,"\n";
		    foreach my $bush_node ( $next_parent->get_all_Descendents ) {
			if( $bush_node->is_Leaf ) {
			    if($bush_node->branch_length < $branch_treshold) {
				if( defined($seen_ids{$bush_node->id}) && $seen_ids{$bush_node->id} == 1 ) {
				    # print "connects to previously shown ",$bush_node->id,"\n";
				    next;
				}
				print $bush_node->id," bl ",$bush_node->branch_length,"\n";
				$seen_ids{$bush_node->id} = 1;
			    }
			}
		    }
		}
	    }
	}
    }
}
