require 5.005;
use ExtUtils::MakeMaker;

my $NAME     = 'ARCT';
my $DISTNAME = "ARCT";
my $VERSION  = "0.9";

check_prerequisites;

WriteMakefile(
	      NAME	  => $NAME,
	      DISTNAME    => $DISTNAME,
	      VERSION     => $VERSION,
	      'dist'      => {  COMPRESS     => 'gzip -9f',
				SUFFIX       => '.gz',
				DIST_DEFAULT => 'all tardist',
			    },
#              'INST_HTMLLIBDIR' => 'doc',
#              'INSTALLHTMLSITELIBDIR' => '$(PREFIX)/doc'

	      );

sub check_prerequisites {
    eval "require Bio::Root::RootI;";
    if ($@) {
		print qq!

		-------------------------------------------------
			            WARNING
		-------------------------------------------------

		 Bioperl does not seem to be installed.

    	 An up to date Bioperl version is required by ARCT.

		 Please install Bioperl BEFORE proceeding with
		 ARCT installation. Go to http://bioperl.org for
		 information on how to obtain and install it.

		-------------------------------------------------

		!;
    }
}

