#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;
use Data::Dumper;

my $xml_check = eval {
	require XML::Simple;
	XML::Simple->import();
	1; 
};
my $db = 'taxonomy';
my $query = 'root[Subtree]';
my $out_file = 'ncbi_tax.entp';
my $out_dir = 'databases';
my $retmax = 100;

my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

my $output = get($url);			#post the esearch URL

#parse WebEnv/Query/Count
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

if (!-d $out_dir) {
	mkdir($out_dir);
}
my $out_path = $out_dir . "/$out_file";
if (-e $out_path) {
	# TODO separate routine to update database
	open(OUT, ">$out_path") || die "Unable to open output file!\n";
	download_tax();
	close OUT;
} else {
	open(OUT, ">$out_path") || die "Unable to open output file!\n";
	download_tax();
	close OUT; 
}

sub download_tax {
	for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
		$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
		$url .= "&retstart=$retstart&retmax=$retmax&rettype=null&retmode=xml";
	    my $efetch_out = get($url);
	    eval{process_fetch($efetch_out)};	#error
	    last;
	}
}

sub update_tax {

}

sub process_fetch {
	my ($data) = @_;
	if ($xml_check) {
		# parse XML data through XML::Simple
		my $xml = new XML::Simple;
		my $parsed_data = $xml->XMLin($data);
		for (my $entry = 0; $entry < $retmax; $entry++) {
			my $sci_name = lc $parsed_data->{Taxon}[$entry]{ScientificName};
			my $tax_num = $parsed_data->{Taxon}[$entry]{TaxId};
			my $lineage = $parsed_data->{Taxon}[$entry]{Lineage};
			print OUT "$sci_name\t$tax_num\t$lineage\n";
		}
	} else {
		# parse XML data through regex if no module present
		my @lineage = ($data =~ /<Lineage>(.+)<\/Lineage>/g);
		my @sci_name = ($data =~ /<ScientificName>(.+)<\/ScientificName>/g);
		my @tax_num = ($data =~ /<TaxId>(\d+)<\/TaxId>/g);
		for (my $ind = 0; $ind < $retmax; $ind++) {
			my $sci = lc $sci_name[$ind];
			print OUT "$sci\t$tax_num[$ind]\t$lineage[$ind]\n";
		}
	}
}