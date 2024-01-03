#
# The ComparativeSystems application
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use File::Temp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run qw(run);
use Cwd;
use Clone;
use P3DataAPI;

my $script = Bio::KBase::AppService::AppScript->new(\&process_compsystems, \&preflight);

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    print STDERR "preflight comparativesystems ", Dumper($params, $app);

    my $token = $app->token();
    my $ws = $app->workspace();

    my $api = P3DataAPI->new();
    my $groups = $params->{genome_groups}; 
    my $numGenomes = 0;
    for my $gg (@$groups) 
    {
        print "$gg\n";
        my $genomes = $api->retrieve_patric_ids_from_genome_group($gg);   
        my $n = @$genomes;
        $numGenomes = $numGenomes + $n;
    }
    my $genomeList = $params->{genome_ids};
    my $glLen = scalar @$genomeList;
    $numGenomes = $numGenomes + $glLen;
    print "$numGenomes genomes\n";

    my $runtime = 0;
    if ($numGenomes < 10) {
        $runtime = 1800;
    } elsif ($numGenomes < 100) {
        $runtime = 3 * 3600;
    } elsif ($numGenomes < 300) {
        $runtime = 6 * 3600;
    } else {
        $runtime = 43200;
    }

    my $pf = {
    cpu => 1,
    memory => '32GB',
    runtime => $runtime,
    storage => 0,
    is_control_task => 0,
    };
    return $pf;
}

sub process_compsystems
{
    my($app, $app_def, $raw_params, $params) = @_;   

    print 'Proc comparative systems ', Dumper($app_def, $raw_params, $params);

    my $token = $app->token();
    my $ws = $app->workspace();

    #
    # Create an output directory under the current dir. App service is meant to invoke
    # the app script in a working directory; we create a folder here to encapsulate
    # the job output.
    #
    # We also create a staging directory for the input files from the workspace.
    # 

    # TODO: may not need a staging directory
    
    # my $cwd = getcwd();
    # my $cwd = File::Temp->newdir( CLEANUP => 1 ); 
    my $cwd = File::Temp->newdir( CLEANUP => 0 ); 
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";

    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";   

    my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    my $dat = { data_api => $data_api };
    my $sstring = encode_json($dat);

    # TODO: is this needed
    my $params_to_app = Clone::clone($params);

    #
    # Write job description.
    #  
    my $jdesc = "$cwd/jobdesc.json";
    open(JDESC, ">", $jdesc) or die "Cannot write $jdesc: $!";
    print JDESC JSON::XS->new->pretty(1)->encode($params_to_app);
    close(JDESC);

    my $parallel = $ENV{P3_ALLOCATED_CPU};

    my @cmd = ("compare_systems","-o",$work_dir,"--jfile", $jdesc);

    warn Dumper (\@cmd, $params_to_app);

    # my $ok = run(\@cmd);
    my $ok = 1;

    if (!$ok)
    {
        die "Command failed: @cmd\n";
    }

    # testing codon tree
    if (1) {
        print "Run codon tree\n";
        my %phylo_fields = (
            'genome_ids' => $params->{genome_ids},
            'genome_group' => $params->{genome_groups}
        );  
        my $output_json = encode_json(\%phylo_fields);
        open(my $file, '>', "$work_dir/file.json") or die "Couldn't open file.json: $!";
        print $file $output_json;
        close($file);

        # my @cmd; 
        # my @options;
        # my $cpu = $ENV{P3_ALLOCATED_CPU};

        # @cmd = ("p3x-build-codon-tree");
        # if ($cpu)
        # {
        #     push(@options, "--threads", $cpu);
        # }

        #push(@options,''    
        
    }

    die "here\n";

    my @output_suffixes = ([qr/\.tsv$/, 'tsv'],[qr/\.json$/, 'json'],[qr/\.txt$/, 'txt']);
    
    my $outfile;
    opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
    # TODO: not sure what this does?
    my @files = sort {$a cmp $b } grep { -f "$work_dir/$_" } readdir(D);

    my $output = 1;
    my $output_dir = "$params->{output_path}/.$params->{output_file}";
    for my $file (@files)
    {
        for my $suf (@output_suffixes)
        {
            if ($file =~ $suf->[0])
            {
                $output = 0;
                my $type = $suf->[1];

                $app->workspace->save_file_to_file("$work_dir/$file", {}, "$output_dir/$file", $type, 1, 
                                                    (-s "$work_dir/$file" > 10_000 ? 1 : 0), #use shock for larger files
                                                    $token);                                                   
            }
        }
    }
}
