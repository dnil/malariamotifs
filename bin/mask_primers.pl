#!/usr/bin/perl -w 

while (<STDIN>) {
    
    # attempt to mask primer artefact aas from alignment

    chomp;

    my $AF = 0;
    $AF += s/ARSFA/arsfa/g;    
    $AF += s/RSFA/rsfa/g;    
    $AF += s/ARGFA/argfa/g;

    $AF > 1 && print STDERR "Found $AF AF.\n";

    my $nDBLf = 0;

    $nDBLf += s/AANYEG/aanyeg/g;
    $nDBLf += s/AANYEA/aanyea/g;
    $nDBLf += s/AAKYEA/aakyea/g;
    $nDBLf += s/AAKYEG/aakyeg/g;
    $nDBLf += s/GSKYEA/gskyea/g;
    $nDBLf += s/GSKYEG/gskyeg/g;
    $nDBLf += s/GSKIEG/gskieg/g;
   
    $nDBLf += s/GTKYEG/gtkyeg/g;

    $nDBLf += s/AKYEG/akyeg/g;
    $nDBLf += s/AIYEG/aiyeg/g;

    $nDBLf += s/SKYEA/skyea/g;
    $nDBLf += s/SKYEG/skyeg/g;
    $nDBLf += s/SNYEG/snyeg/g;
    
    $nDBLf > 1 && print STDERR "Found $nDBLf nDBLf.\n";

    my $nDBLr = 0;
    

    $nDBLr += s/REDWW/redww/g; # 5-7
    $nDBLr > 1 && print STDERR "Found $nDBLr nDBLr.\n";

    my $BR = 0;
    
    $BR += s/WFDEW/wfdew/g;
    $BR += s/WFDEM/wfdem/g;
    $BR += s/WFDEL/wfdel/g;
    $BR += s/WFEEW/wfeew/g;
    $BR += s/WFEGM/wfegm/g;

    $BR > 1 && print STDERR "Found $BR BR.\n";

    print $_, "\n";


}
