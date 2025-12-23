#!/usr/bin/perl
use strict;
use warnings;

# -------------------------
# Usage:
#   perl softclip_assembly_blat_sv.pl sample.bam chr12 11648801 11836600
# -------------------------

# Configuration / inputs
my $bam_file = $ARGV[0] or die "Usage: $0 <bam> <chrom> <region_start> <region_end>\n";
my $chrom    = $ARGV[1] or die "Usage: $0 <bam> <chrom> <region_start> <region_end>\n";
my $region_start = $ARGV[2] or die "Usage: $0 <bam> <chrom> <region_start> <region_end>\n";
my $region_end   = $ARGV[3] or die "Usage: $0 <bam> <chrom> <region_start> <region_end>\n";

my $upstream_distance   = 100;  # adjust as needed
my $downstream_distance = 100;  # adjust as needed

# Read-quality filters
my $MIN_MAPQ        = 20;   # mapping quality cutoff
my $MIN_BQ          = 20;   # base quality cutoff (Phred)
my $MIN_CLIP        = 10;   # minimum soft-clip length
my $MAX_LOW_BQ_FRAC = 0.2;  # allow up to 20% low-BQ bases within the clip

# Assembly settings
my $cap3_path = "cap3";             # CAP3 executable (or full path)
my $MIN_READS_FOR_ASSEMBLY = 5;     # skip assembly if too few clipped reads

# BLAT / contig alignment coverage threshold
my $MIN_CONTIG_ALN_FRAC = 0.98;     # at least 98% of contig must be aligned overall

# Reference / tools
my @samples = split(/\//, $bam_file);
my $sample  = pop(@samples);
$sample =~ s/\.bam$//;

# EDIT this path as needed:
my $reference_fasta = "/research/groups/geelegrp/projects/PedDep_DataScience/common/WGS/PedDepAcc/bwang/genome/hg38.analysisSet.fa";
my $blat_path       = "blat";

my $temp_dir = "temp_blat_$sample.$ARGV[1].$ARGV[2].$ARGV[3]";
mkdir $temp_dir unless -d $temp_dir;

# Define regions (4 windows around two boundaries)
my %regions = (
    "boundary1_upstream"   => [$region_start - $upstream_distance,         $region_start - 1],
    "boundary1_downstream" => [$region_start,                               $region_start + $downstream_distance - 1],
    "boundary2_upstream"   => [$region_end   - $upstream_distance,         $region_end   - 1],
    "boundary2_downstream" => [$region_end,                                 $region_end   + $downstream_distance - 1],
);

# -------------------------
# Extract reads with high-quality soft clips (high-quality only)
# Returns arrayref of FULL read sequences (one entry per read)
# -------------------------
sub extract_soft_clipped_reads {
    my ($region_name, $chrom, $start, $end) = @_;

    # Avoid invalid coords
    $start = 1 if $start < 1;
    $end   = 1 if $end   < 1;

    my $region_bam = "$temp_dir/$region_name.bam";
    my $region_sam = "$temp_dir/$region_name.sam";

    # Extract reads in region
    system("samtools view -b $bam_file $chrom:$start-$end > $region_bam") == 0
        or die "samtools view failed for $region_bam";

    # Convert BAM to SAM
    system("samtools view $region_bam > $region_sam") == 0
        or die "samtools view failed for SAM conversion: $region_sam";

    open(my $fh, "<", $region_sam) or die "Cannot open $region_sam: $!";

    my @sequences;
    my %seen_read;

    while (my $line = <$fh>) {
        next if $line =~ /^@/;
        chomp $line;

        my @fields = split("\t", $line);
        next if @fields < 11;

        my ($flag, $mapq, $cigar, $seq, $qual) = @fields[1, 4, 5, 9, 10];

        next if !defined $cigar || $cigar eq "*";
        next if !defined $seq   || $seq   eq "*";
        next if !defined $qual  || $qual  eq "*";

        # Skip unmapped / secondary / supplementary
        next if ($flag & 0x4);
        next if ($flag & 0x100);
        next if ($flag & 0x800);

        # MAPQ filter
        next if $mapq < $MIN_MAPQ;

        # A read can have left and/or right soft clips
        my @clips;

        # Identify left soft-clip
        if ($cigar =~ /^(\d+)S/) {
            my $clip_len = $1;
            if ($clip_len >= $MIN_CLIP) {
                my $clip_seq  = substr($seq,  0, $clip_len);
                my $clip_qual = substr($qual, 0, $clip_len);
                push @clips, [$clip_seq, $clip_qual] if length($clip_seq) == length($clip_qual);
            }
        }

        # Identify right soft-clip
        if ($cigar =~ /(\d+)S$/) {
            my $clip_len = $1;
            if ($clip_len >= $MIN_CLIP) {
                my $clip_seq  = substr($seq,  -$clip_len);
                my $clip_qual = substr($qual, -$clip_len);
                push @clips, [$clip_seq, $clip_qual] if length($clip_seq) == length($clip_qual);
            }
        }

        # Evaluate base-quality in clipped bases; if any soft-clip
        # passes quality filters, keep the FULL READ (not only the clip)
        my $has_good_softclip = 0;
        my $qname = $fields[0];

        for my $c (@clips) {
            my ($clip_seq, $clip_qual) = @$c;
            next if length($clip_qual) == 0;

            my $low_bq = 0;
            foreach my $q (split //, $clip_qual) {
                my $phred = ord($q) - 33;
                $low_bq++ if $phred < $MIN_BQ;
            }

            next if ($low_bq / length($clip_qual)) > $MAX_LOW_BQ_FRAC;

            # This soft-clipped side is good enough; mark the read
            $has_good_softclip = 1;
            last;
        }

        # If this read has at least one good soft-clipped side,
        # store the FULL read sequence once per region
        if ($has_good_softclip) {
            next if $seen_read{$qname}++;
            push @sequences, $seq;
        }
    }

    close($fh);
    return \@sequences;
}

# -------------------------
# Write sequences to FASTA
# -------------------------
sub write_fasta {
    my ($region, $seqs_ref) = @_;
    my $fasta_file = "$temp_dir/$region.fasta";

    open(my $out, ">", $fasta_file) or die "Cannot write to $fasta_file: $!";
    my $idx = 1;
    for my $seq (@$seqs_ref) {
        print $out ">read_${region}_$idx\n$seq\n";
        $idx++;
    }
    close($out);

    return $fasta_file;
}

# -------------------------
# Assemble with CAP3
# Returns path to FASTA to BLAT: contigs if present, else singlets, else original
# -------------------------
sub assemble_with_cap3 {
    my ($region, $input_fasta) = @_;

    my $log = "$temp_dir/$region.cap3.log";
    my $cmd = "$cap3_path $input_fasta > $log 2>&1";
    system($cmd) == 0 or die "CAP3 failed for $region. Check $log\n";

    # CAP3 produces: <input>.cap.contigs and <input>.cap.singlets
    my $contigs  = "$input_fasta.cap.contigs";
    my $singlets = "$input_fasta.cap.singlets";

    if (-s $contigs) {
        return $contigs;
    } elsif (-s $singlets) {
        return $singlets;
    } else {
        return $input_fasta;
    }
}

# =======================================================
# SV classification helpers
# =======================================================

# classify two genomic anchors into SV type and breakpoint strings
sub classify_sv_pair {
    my ($hit1, $hit2) = @_;

    my $svtype;
    my ($chr1, $chr2) = ($hit1->{tName}, $hit2->{tName});

    my $bp1 = join(".", $hit1->{tName}, $hit1->{tEnd},   $hit1->{strand});
    my $bp2 = join(".", $hit2->{tName}, $hit2->{tStart}, $hit2->{strand});

    if ($chr1 ne $chr2) {
        $svtype = "CTX";  # interchromosomal translocation
        return ($svtype, $bp1, $bp2);
    }

    # same chromosome
    my $pos1 = $hit1->{tEnd};
    my $pos2 = $hit2->{tStart};

    my $qgap = abs($hit1->{qEnd} - $hit2->{qStart});
    my $tdiff = abs($pos2 - $pos1);
    my $buffer = 10;

    if ($hit1->{strand} ne $hit2->{strand}) {
        # orientation change on same chromosome: inversion breakpoint
        $svtype = "inversion";
        return ($svtype, $bp1, $bp2);
    }

    if ($tdiff > $qgap + $buffer) {
        # genome gap >> contig gap -> deletion in sample (sequence missing)
        $svtype = "deletion";
    } elsif ($qgap > $tdiff + $buffer) {
        # contig gap >> genome gap -> insertion in sample (extra sequence)
        $svtype = "insertion";
    } else {
        $svtype = "indel_or_complex";
    }

    return ($svtype, $bp1, $bp2);
}

# -------------------------
# Helper: compute union coverage of query intervals across hits
#         returns total covered length on the contig
# -------------------------
sub compute_query_coverage {
    my ($hits_ref) = @_;
    return 0 unless @$hits_ref;

    my @intervals = map { [ $_->{qStart}, $_->{qEnd} ] } @$hits_ref;
    @intervals = sort { $a->[0] <=> $b->[0] } @intervals;

    my $cov = 0;
    my ($cur_s, $cur_e) = @{$intervals[0]};

    for my $i (1 .. $#intervals) {
        my ($s, $e) = @{$intervals[$i]};
        if ($s <= $cur_e) {
            # overlap or adjacency: extend current interval
            $cur_e = $e if $e > $cur_e;
        } else {
            # disjoint interval
            $cov += $cur_e - $cur_s;
            ($cur_s, $cur_e) = ($s, $e);
        }
    }

    $cov += $cur_e - $cur_s if @intervals;

    return $cov;
}

# -------------------------
# Run BLAT, then detect SVs from all PSL lines per contig
# Only consider contigs where >=98% of qSize is aligned
# -------------------------
sub blat_and_classify_sv {
    my ($region, $query_fasta) = @_;

    my $psl_output = "$temp_dir/$region.psl";

    my $cmd = "$blat_path -fine $reference_fasta $query_fasta $psl_output";
    system($cmd) == 0 or warn "BLAT failed for $region\n";

    return unless -s $psl_output;

    open(my $fh, "<", $psl_output) or die "Cannot open $psl_output: $!";

    my %hits_by_qname;

    while (my $line = <$fh>) {
        next if $line =~ /^\s*#/;
        next if $line =~ /^\s*$/;

        chomp $line;
        my @f = split(/\t/, $line);

        next unless @f >= 21;

        my $hit = {
            matches    => $f[0] + 0,
            misMatches => $f[1] + 0,
            repMatches => $f[2] + 0,
            nCount     => $f[3] + 0,
            qNumInsert => $f[4] + 0,
            tNumInsert => $f[6] + 0,
            strand     => $f[8],
            qName      => $f[9],
            qSize      => $f[10] + 0,
            qStart     => $f[11] + 0,
            qEnd       => $f[12] + 0,
            tName      => $f[13],
            tStart     => $f[15] + 0,
            tEnd       => $f[16] + 0,
            raw        => $line,
        };

        # Only keep reasonably good hits
        next if $hit->{matches} < 20;

        push @{$hits_by_qname{ $hit->{qName} }}, $hit;
    }
    close($fh);

    # For each contig (qName)
    for my $qname (sort keys %hits_by_qname) {
        my @hits = sort { $a->{qStart} <=> $b->{qStart} } @{$hits_by_qname{$qname}};
        next if @hits < 2;

        # Require that overall aligned length on contig >= 98% of qSize
        my $qsize = $hits[0]->{qSize} || 0;
        next if $qsize <= 0;

        my $cov = compute_query_coverage(\@hits);
        my $frac = $cov / $qsize;

        # Skip SV calls for contigs that are not well-covered
        next if $frac < $MIN_CONTIG_ALN_FRAC;

        for (my $i = 0; $i < @hits - 1; $i++) {
            my $h1 = $hits[$i];
            my $h2 = $hits[$i+1];

            my ($svtype, $bp1, $bp2) = classify_sv_pair($h1, $h2);

            print join("\t",
                $region,
                $qname,
                $svtype,
                $bp1,
                $bp2,
                $h1->{raw},
                $h2->{raw}
            ), "\n";
        }
    }
}

# =======================================================
# Main: for each region, extract soft-clipped FULL reads,
#       assemble with CAP3, then BLAT and classify SVs.
# =======================================================

my %region_query_fasta;

for my $region (keys %regions) {
    my ($start, $end) = @{$regions{$region}};

    my $seqs_ref = extract_soft_clipped_reads($region, $chrom, $start, $end);
    my $count    = scalar(@$seqs_ref);

    next unless $count;

    print "Region $region has $count softclip-supporting reads\n";

    my $raw_fasta = write_fasta($region, $seqs_ref);

    my $query_fasta = $raw_fasta;
    if ($count >= $MIN_READS_FOR_ASSEMBLY) {
        $query_fasta = assemble_with_cap3($region, $raw_fasta);
        print "Region $region assembly output used for BLAT: $query_fasta\n";
    } else {
        print "Region $region has too few sequences for assembly; BLAT raw reads.\n";
    }

    $region_query_fasta{$region} = $query_fasta;
}

for my $region (sort keys %region_query_fasta) {
    blat_and_classify_sv($region, $region_query_fasta{$region});
}

# Cleanup (optional)
# system("rm -rf $temp_dir");

