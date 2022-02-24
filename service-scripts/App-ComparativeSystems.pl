#
# The FastQ Utils application.
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::ReadSet;

use strict;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run qw(run);
use Cwd;
use Clone;

my $script = Bio::KBase::AppService::AppScript->new(\&process_fastq, \&preflight);

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    print STDERR "preflight fastqc ", Dumper($params, $app);

    my $token = $app->token();
    my $ws = $app->workspace();

    my $readset;
    eval {
	$readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params);
    };
    if ($@)
    {
	die "Error parsing assembly parameters: $@";
    }

    my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($ws);

    if (!$ok)
    {
	die "Reads as defined in parameters failed to validate. Errors:\n\t" . join("\n\t", @$errs);
    }
    print STDERR "comp=$comp_size uncomp=$uncomp_size\n";

    my $est_uncomp = $comp_size / 0.75 + $uncomp_size;

    my $est_time = int($est_uncomp * 1e-6 * 3.0);

    my $est_cpu = 8;
    my $est_ram = '32G';

    if ($est_time < 3600)
    {
	$est_time = 3600;
    }
    elsif ($est_time > 3600 * 2)
    {
	$est_ram = '128G';
    }

    return {
	cpu => $est_cpu,
	memory => $est_ram,
	runtime => $est_time,
    };
}


sub process_fastq
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Proc fastq utils ", Dumper($app_def, $raw_params, $params);

    my $token = $app->token();
    my $output_folder = $app->result_folder();

    #
    # Create an output directory under the current dir. App service is meant to invoke
    # the app script in a working directory; we create a folder here to encapsulate
    # the job output.
    #
    # We also create a staging directory for the input files from the workspace.
    #

    my $cwd = getcwd();
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";

    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

    my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    my $dat = { data_api => $data_api };
    my $sstring = encode_json($dat);

    #
    # Read parameters and discover input files that need to be staged.
    #
    # Make a clone so we can maintain a list of refs to the paths to be
    # rewritten.
    #
    my %in_files;

    my $params_to_app = Clone::clone($params);
    my @to_stage;

    for my $read_tuple (@{$params_to_app->{paired_end_libs}})
    {
	for my $read_name (keys %{$read_tuple})
	{
	   if($read_name == "read1" || $read_name == "read2")
           {
	       my $nameref = \$read_tuple->{$read_name};
	       $in_files{$$nameref} = $nameref;
	       push(@to_stage, $$nameref);
           }
        }
    }
    for my $read_tuple (@{$params_to_app->{single_end_libs}})
    {
	for my $read_name (keys %{$read_tuple})
	{
	   if($read_name == "read")
           {
	       my $nameref = \$read_tuple->{$read_name};
	       $in_files{$$nameref} = $nameref;
	       push(@to_stage, $$nameref);
           }
        }
    }

    my $staged = {};
    if (@to_stage)
    {
	warn Dumper(\%in_files, \@to_stage);
	$staged = $app->stage_in(\@to_stage, $stage_dir, 1);
	while (my($orig, $staged_file) = each %$staged)
	{
	    my $path_ref = $in_files{$orig};
	    $$path_ref = $staged_file;
	}
    }

    #
    # Write job description.
    #
    my $jdesc = "$cwd/jobdesc.json";
    open(JDESC, ">", $jdesc) or die "Cannot write $jdesc: $!";
    print JDESC JSON::XS->new->pretty(1)->encode($params_to_app);
    close(JDESC);

    my $parallel = $ENV{P3_ALLOCATED_CPU};
    my $override = {
	fastqc => { -p => $parallel},
	trim_galore => {-p => $parallel},
	bowtie2 => {-p => $parallel},
	hisat2 => {-p => $parallel},
	samtools_view => {-p => $parallel},
	samtools_index => {-p => $parallel},
    samtools_sort => {-p => $parallel}
    };

    my @cmd = ("p3-fqutils",
	       "--jfile", $jdesc,
	       "--sstring", $sstring,
	       "-p", encode_json($override),
	       "-o", $work_dir);

    warn Dumper(\@cmd, $params_to_app);

    my $ok = run(\@cmd);
    # my $ok = run(\@cmd,
	# 	 ">", "$work_dir/fqutils.out.txt",
	# 	 "2>", "$work_dir/fqutils.err.txt");
    if (!$ok)
    {
        # opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
        # $app->workspace->save_file_to_file("$work_dir/fqutils.out.txt", {}, "$output_folder/fqutils.out.txt", "txt", 1,
        #                     (-s "$work_dir/fqutils.out.txt" > 10_000 ? 1 : 0), # use shock for larger files
        #                     $token);
        # $app->workspace->save_file_to_file("$work_dir/fqutils.err.txt", {}, "$output_folder/fqutils.err.txt", "txt", 1,
        # (-s "$work_dir/fqutils.err.txt" > 10_000 ? 1 : 0), # use shock for larger files
        # $token);
	    die "Command failed: @cmd\n";
    }

    my @output_suffixes = ([qr/\.bam$/, "bam"],
			   [qr/\.fq\.gz$/, "reads"],
			   [qr/\.fq\.1.gz$/, "reads"],
			   [qr/\.fq\.2.gz$/, "reads"],
			   [qr/\.bai$/, "bai"],
			   [qr/\.html$/, "html"],
			   [qr/\.fastq\.gz$/, "reads"],
			   [qr/\.txt$/, "txt"]);

    my $outfile;
    opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
    my @files = sort { $a cmp $b } grep { -f "$work_dir/$_" } readdir(D);

    my $output=1;
    for my $file (@files)
    {
	for my $suf (@output_suffixes)
	{
	    if ($file =~ $suf->[0])
	    {
 	    	$output=0;
		my $path = "$output_folder/$file";
		my $type = $suf->[1];

		$app->workspace->save_file_to_file("$work_dir/$file", {}, "$output_folder/$file", $type, 1,
					       (-s "$work_dir/$file" > 10_000 ? 1 : 0), # use shock for larger files
					       $token);
	    }
	}
    }

    #
    # Clean up staged input files.
    #
    while (my($orig, $staged_file) = each %$staged)
    {
	unlink($staged_file) or warn "Unable to unlink $staged_file: $!";
    }

    return $output;
}
