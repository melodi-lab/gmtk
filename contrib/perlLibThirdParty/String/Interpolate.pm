use strict;
use warnings;

package String::Interpolate;
our $VERSION = 0.3;
use Carp qw( croak );

=head1 NAME

String::Interpolate - Wrapper for builtin the Perl interpolation engine.

=head1 SYNOPSIS
    
    # Functional interface
    use String::Interpolate qw( safe_interpolate interpolate ); 
    our($GREET) = 'Hello'; # Cannot be lexical 
    print interpolate( '$GREET $1\n', [ 'world' ] );

    # Object interface
    use String::Interpolate;
    my $who;
    my $template = new String::Interpolate { WHO => \$who };
    $template->{TIME} = sub () { localtime }; # Tie $TIME to localtime()
    $template->( [ qw( now it ) ] ); # Set $1, $2
    $template->[3] = 'is'; # Sets $3
    $who = 'old friend';
    $template->( '$REV{olleH} $WHO, $2 $3 $1 $TIME$_' ); # Set string to process
    $template->{REV} = sub { reverse @_ }; # Tie %REV to reverse()
    $_ = '.';
    print "$template\n"; # Perform interpolation

    # Peform the interpolation in a Safe compartment.
    my $replace = safe String::Interpolate '\u\L$1'; 
    my $search = qr/(\w+)/;
    $_ = "HELLO world\n";
    s/$search/$replace/eg; # /e supresses optimisation
    print;

=head1 DESCRIPTION

C<String::Interpolate> provides a neat interface to the solution to
that perenial Perl problem - how to invoke the Perl string
interpolation engine on a string contained in a scalar variable.

A C<String::Interpolate> object encapsulates a string and a context in
which it should be subjected to Perl interpolation.  In the
simplest, default, case the context is simply the namespace (package)
from which the constructor was called.

A C<String::Interpolate> object may hold a reference to an array and
hashes that will be used to populate the special variables $1 etc and
some package variables respectively prior to each interpolation.  

In general special globally global variables such as $_ can be used in
the interpolation, the exception being @_ which is always empty during
the interpolation.

The interpolated string is processed with strictures and warnings
enabled excluding 'strict vars' and 'warnings uninitialized' so that
interpolating undefined variables will be silently ignored.  This
behaviour can be altered using the pragma() method.

Because the Perl string interpolation engine can call arbitrary Perl
code you do not want to want to use it on strings from untrusted
sources without some precautions.  For this reason
C<String::Interpolate> objects can be made to use C<Safe>
compartments.  This is, of course, only as safe as Safe and you are
advised to read "WARNING" section of the Safe documentation.

When interpolating in a Safe compartment package symbols are imported
using tied wrapper variables so that their values cannot be
interpreted as references and such that they cannot be used to alter
the values outside the compartment.  This behaviour can be suppressed
by the unsafe_symbols() method.  Note that if you want to import tied
variable or variables containing references to objects that use
overloading into a Safe compartment then you will need to do a lot of
fancy footwork unless you use safe_hole() method.

By default *_ is shared by Safe compartments and could potentially
allow the compartment to leak.  The $_ and %_ variables are therefore
subjected to the same similar precautions to imported symbols.  This
behaviour can be suppressed using the unsafe_underscore() method.

Perl string interpolation can, of course, throw exceptions.  By
default String::Interpolate objects do not catch (or rethrow) these
exceptions when working in a simple namespace and do trap them when
working in a Safe compartment.  This behaviour can be overriden by the
trap() or pragma() methods.  If an exception during interpolation is
trapped then undef will be returned as the result of the
interpolation and $@ will hold the exception in the usual way.

When taint checking enabled, attempting to perform interpolation
(using eval()) on a tainted string would naturally fail.  However,
when using a Safe compartment, String::Interpolate will strip the
tainting off of the string prior to interpolation and put it back
afterwards.  Also String::Interpolate will taint any arguments
passed to callback functions called as the result of performing
interpolation on a tainted string.  Note that due to the mechanism
used to assign $1 et al they can never be tained even if the values in
the array being used to set them are tainted.

By default C<String::Interpolate> does not export any subroutines but
as a concession to programmers who prefer not to explicitly use
objects the functions interpolate() and safe_interpolate() are
exportable.

=cut

# Must appear before any file-scoped lexicals
sub reval { no strict 'vars'; eval $_[0] }

sub prevent_blessed_error_hack () {
    return unless ref $@;
    no strict 'refs';
    no warnings 'redefine';
    local *{"@{[ref $@]}::DESTROY"} = sub {};
    $@ = 'Blessed error from Safe compartment';
}

# During Carp::confess stack dumps we don't want to exec()
# %dbgpkg is a package variable as callers may want to manipulate it.

our %dbgpkg = (
	     Carp => 1,
	     );

our $taint_flag = '';
our $safe_hole;

my %type_from_prefix = (
			"\$" => 'SCALAR',
			'@' => 'ARRAY',
			'%' => 'HASH',
			);

use overload 
    '""' => sub { 
	my $self = shift;
	$dbgpkg{caller()} ? overload::StrVal($self) : $self->exec;
    },
    'cmp' => sub { my ($l,$r) = @_; $l->exec cmp $r },
    '@{}' => sub { tie my @a, 'String::Interpolate::AsArray', @_; \@a },
    '%{}' => 'ashash',
    '&{}' => sub { my $self=shift; sub { $self->exec(@_) } };


use base 'Exporter';
our(@EXPORT_OK) = qw ( interpolate safe_interpolate );
my $pkgcount;

=head2 Principle methods

=over 4

=item new

Simple constructor.  Creates a empty String::Interpolate object bound
to the caller's namespace and then modifies the object by passing any
arguments to the exec() method.  Returns a the object.

If called as an instance method new() clones the object.  Be aware,
however, that this is a shallow cloning and if array or hash reference
arguments have been passed to the object the parent and clone will
continue to use the same array or hashes until one or other is passed
a new argument.

Most of the other methods in String::Interpolate will implicitly call
new() if called as class methods.

=cut

my %preset_pragma = (
   NOWARN => 'unimport warnings qw(uninitialized)',
   WARN => '',
   FATAL => 'import warnings FATAL => qw(uninitialized); import strict qw(vars)',
 );
		      
sub new {
    my $class = shift;
    my $self;
    if ( ref $class ) {
	# Clone
	$self = bless \ { %$$class }, ref $class;
	delete @$$self{'tmppkg','pkg','code'} if $$self->{tmppkg};
	delete $$self->{safe} if $$self->{implicit_safe};
    } else {
	my $calldepth = 0;
	my $defpgk;
	do { $defpgk = caller($calldepth++) } while $defpgk->isa( __PACKAGE__ );
        $self = bless \ {
	    defpgk => $defpgk,
	    pkg => $defpgk,
	    pragmas => $preset_pragma{NOWARN},
	}, $class;
    }
    $self->exec(@_);
    $self;
}

=item safe

Alternative constuctor to create a String::Interpolate object that
uses an automatically allocated temporary Safe compartment.  The
automatically allocated Safe compartment will have the default opcode
mask but with the 'bless' opcode denied as this can be used to execute
code outside the compartment by putting it in DESTROY methods.  The
'tie' opcode is also denied although I'm not sure if it really can be
exploited in this way.  There is no point explicitly passing a package
or existing safe compartment to this constructor as it will be ignored.

The argument list is passed to exec() as in new().

The safe() method can also be called on an existing object in which
case it instructs the object to forget its current Safe compartment or
namespace and use an automatically allocated temporary Safe
compartment henceforth.

=cut

sub safe {
    my $self = shift;
    $self = $self->new(@_) unless ref $self;
    $self->free_tmppkg;
    delete @$$self{'pkg','explicit_pkg','safe'};
    $$self->{implicit_safe}++;
    require Safe;
    $self;
}

=item exec

This it the guts of the implementation but it it rarely needs to be
called explicitly as it can be more elegantly called implicitly by
using the String::Interpolate object in a string or CODE reference
context.  The following are equivalent pairs:

    my $interpolated_string = $interpolate_object->exec;
    my $interpolated_string = "$interpolate_object";

    my $interpolated_string = $interpolate_object->exec(LIST);
    my $interpolated_string = $interpolate_object->(LIST);

The exec() method modifies the object according the argument list.
Then, if called in a non-void context, returns the result of the
interpolation.  Note that the modifications are persistent.  This
persistence can be avoided by creating a transient clone using the
new() method.

    my $string = $inter->(LIST);      # $inter changed
    my $string = $inter->new->(LIST); # $inter unchanged

Also, if exec() is called as a class method then it acts on a
temporary String::Interpolate object which is immediately destroyed.

The elements of the argument list are interpreted according to their
type as listed below.  If this mechanism does not provide sufficient
flexibility in manipulating the symbol table you can, of course,
manipulate it directly too.

=over 4

=item ARRAY reference

Tells the object to use this array to populate the special variables
$1 and so on.  The object holds a reference to the array itself and
will use the values that are in the array at the time of
interpolation.  This ARRAY reference is exposed via the positionals()
method.  The array can also be modified by using the
String::Interpolate object in an ARRAY reference context.  Note,
however, that the String::Interpolate object used in an ARRAY
reference context does not refer to the array itself but to a
STORE-only tied array whose subscripts are offset by one such that
$interpolate_object->[1] corresponds to
$interpolate_object->positionals->[0] and hence the value that will be
interpolated for $1.

=item HASH reference

Tells the object to use this hash to populate some package variables
immediately prior to each interpolation.  The object holds a reference
to the hash itself and will use the values that are in the hash at the
time of interpolation. 

After the object has been instructed to populate package variables in
this way it will no longer default to using the namespace from which
the constructor was called and will instead auto-allocate a temporary
one unless told to do otherwise.

If multiple hash reference arguments are specified in a single call to
exec() then each hash in turn will be processed prior to each
interpolation.  However, whenever a exec() is passed one or more hash
references it forgets any previous hashes and deletes any
auto-allocated temporary package or safe compartment.

The keys of the hash should be unqualified Perl identifiers that will
determine the entries in the package symbol to be modified.  Which slot
in the symbol table entry is modified is determined by the values'
types as follows:

=over 4

=item ARRAY reference

Set the symbol table entry's ARRAY slot.

=item HASH reference

Set the symbol table entry's HASH slot.

=item SCALAR reference

Set the symbol table entry's SCALAR slot.

=item CODE reference with prototype ()

Set the symbol table entry's SCALAR slot to point to an new tied
scalar with a FETCH method that calls the referenced code.  

Note that if interpolation is taking place inside a Safe compartment
the callback will, by default, simply be called from within the
compartment.  The callback code will execute with a false symbol table
root so it will not be able to use any packages from the real symbol
table root.  This limitation can be overcome by using the safe_hole()
method.

=item CODE reference with prototype ($) or no prototype

Set the symbol table entry's HASH slot to point to an new tied
hash with a FETCH method that calls the referenced code.

See above for limitations if the callback is called from interpolation
taking place in a Safe compartment.

The argument passed to the callback will be stringified.  It may seem
like a nice idea to be able to pass multiple arguments using an ARRAY
reference but unfortunately this could open up security problems when
passing arguments out of a Safe compartment via a Safe::Hole.

=item Anything else 

Set the symbol table entry's SCALAR slot to point
scalar containing the value.

=back

Note that since the String::Interpolate object stores a reference to
the hash and updates the symbol table prior to each interpolation,
changes in the hash will be reflected in subsequent interpolations.
However, if items in the hash are deleted or changed to a different
type then the previously created symbol table entries may persist.
This can be overcome by calling the safe() or package() methods.

To simplify modifying the hash, a String::Interpolated object used in
a HASH reference context will return a reference to the last hash
argument passed to object, implicitly calling exec({}) first if
necessary.

    my %h = ( A => 1 );
    my $i = new String::Interpolate \%h;
    $i->{B} = 2;  # $h{B} = 2

=item GLOB or GLOB reference

Instruct the object to perform interpolation in the namespace defined
by the GLOB.  For example the argument *Q:: would mean that the string
should be interpolated in the context of the package Q.  The trailing
'::' may be omitted.  

Passing a package argument to the object causes it to stop using a
Safe compartment if it previously was doing so.  If you want safe
execution in a specific namespace then you need to explicitly constuct
Safe object bound to the given namespace and pass that.

Once a String::Interpolate object has been explicitly bound to a
namespace it will continue to use that namespace even if the
String::Interpolate object has been (or is subsequently) passed a hash
reference argument.  In this case the symbols will be created/updated
in the namespace prior to each interpolation and will persist
afterwards.

See also the package() method.

=item Safe object

Instruct the object to perform interpolation in the given Safe
compartment.  Passing a Safe object argument to the
String::Interpolate object causes it to stop using a specified
namespace if it previously was doing so.  If you choose to pass an
explicit Safe object you should deny the 'bless' and 'tie' opcodes for
the reasons discussed under the safe() method.

Once a String::Interpolate object has been explicitly bound to a Safe
object it will continue to use that object even if the
String::Interpolate object has been (or is subsequently) passed a hash
reference argument.  In this case the symbols will be created/updated
in the namespace associated with the Safe object prior to each
interpolation and will persist afterwards.

See also the safe() method.

=item Safe::Hole object

Equivalent to calling the safe_hole() method with the same argument.

=item SCALAR reference

The referenced scalar is passed to the pragma() method.

=item Anything else

Use the stringified value of the argument as the string on which to
perform interpolation.

=back

=cut

sub exec {
    my $self = shift;
    $self = $self->new unless ref $self;
    my $seenmap;

    for ( @_ ) {
	if ( ref eq 'ARRAY' ) {
	    $$self->{pos} = $_;
	} elsif ( ref eq 'HASH' ) {
	    my $map = \$$self->{map};
	    if ( !$seenmap++ && $$map && @$$map ){
		$$map = [];
		$self->free_tmppkg;
	    }
	    push @$$map => $_;
	} elsif ( ref $_ eq 'SCALAR' ) {
	    $self->pragma($$_);
	} elsif ( ref $_ eq 'GLOB' || ref \$_ eq 'GLOB' ) {
	    $self->package($_);
	} elsif ( ref && $_->isa('Safe::Hole') ) {
	    $$self->{safe_hole} = $_;
	} elsif ( ref && $_->isa('Safe') ) {
	    $self->free_tmppkg;
	    delete $$self->{pkg};
	    delete $$self->{implicit_safe};
	    delete $$self->{lexicals};
	    $$self->{safe} = $_;
	    $$self->{trap} = 1 unless defined $$self->{trap};
	} else {
	    $$self->{string} = "$_";
	    delete $$self->{code};
	}
    }
    return unless defined wantarray;

    @_ = ();
    local $_ = $_;

    my $string = $$self->{string};
    my $pos = $$self->{pos};
    my $pkg = $$self->{pkg};
    my $safe = $$self->{safe};
    my $code = $$self->{code};

    if ( $$self->{implicit_safe} && !$safe ) {
	$safe = $$self->{safe} = Safe->new;
	$safe->deny('tie','bless');
    }

    my $dlm = '_aaa';

    if ( defined $string && !$code || $pos ) {
	my $cat = join '' => $string, @{ $pos || [] };
	$dlm++ while -1 < index $cat, $dlm;
    }

    ( join $dlm => @$pos ) =~ /^@{[ join $dlm => ('(.*)') x @$pos ]}$/ 
	or die 'Unexpected pattern match failure initialising $1 et al'
	    if $pos;
 
    if ( $pkg && $pkg eq 'Safe') {
	require Safe;
	$safe = Safe->new;
    }

    $pkg = $safe->root if $safe;

    local $_ = do { no warnings 'uninitialized'; "$_"},
    local *_ = %_ ? String::Interpolate::Func->wrap_hash('_',\%_) : {}
    if $safe && ! $$self->{unsafe_underscore};

    my $safe_symbols = $safe && ! $$self->{unsafe_symbols};

    # use PadWalker qw( peek_my ); use Data::Dumper; die Dumper peek_my(2);
    
    my @pad_map;

    if ( $$self->{lexicals} ) {
	my $depth = 1;
	$depth++ while caller($depth)->isa(__PACKAGE__);
	# die "$depth ". scalar(caller($depth));
	require PadWalker;
	my $pad = PadWalker::peek_my($depth+1);
	# use Data::Dumper; die Dumper $pad;
	while ( my ( $k,$v ) = each %$pad ) {
	    $k =~ s/^([@%\$])//
		or die "$k does not start with \$, \@ or \%";
	    $v = *$v{$type_from_prefix{$1}} if ref $v eq 'GLOB';
	    push @pad_map => { $k => $v };
	}
    }

    for ( @pad_map, @{$$self->{map}} ) {
	$pkg ||= $$self->{tmppkg} ||= __PACKAGE__ . '::' . ++$pkgcount;
	while ( my ( $k,$v ) = each %$_ ) {
	    no strict 'refs';
	    *{"${pkg}::$k"} = do {
		if ( ref $v eq 'HASH' ) {
		    if ( $safe_symbols ) {
		        String::Interpolate::Func->wrap_hash($k,$v);
		    } else {
			$v;
		    }
		} elsif ( ref $v eq 'CODE' ) {
		    my $p = prototype($v);
		    if ( defined $p && ! $p ) {
			my $unimplemented = sub {
			    croak "\$$k tied scalar is FETCH-only within String::Interpolate";
			};
			tie my $s, 'String::Interpolate::Func', {
			    FETCH => $v,
			    STORE => $unimplemented,
			};
			\$s;
		    } elsif ( $p && $p ne "\$" ) {
  		        croak "Invalid prototype ($p) for interpolated function $k";
		    } else {
			my $unimplemented = sub {
			    die "%$k tied hash is FETCH-only within String::Interpolate";
			};
			tie my %h, 'String::Interpolate::Func', {
			    FETCH => $v,
			    STORE => $unimplemented,
			    DELETE => $unimplemented,
			    FIRSTKEY => $unimplemented,
			    NEXTKEY => $unimplemented,
			};
			\%h;
		    }
		} elsif ( ref $v eq 'ARRAY' ) {
		    if ( $safe_symbols ) {
			my $unimplemented = sub {
			    die "\@$k is read-only within String::Interpolate";
			};
			tie my @a, 'String::Interpolate::Func', {
			    FETCH => sub { "$v->[$_[0]]" },
			    STORE => $unimplemented,
			    DELETE => $unimplemented,
			    FETCHSIZE => sub { scalar @$v },
			};
			\@a;
		    } else {
			$v;
		    }
		} elsif ( ref $v eq 'SCALAR' ) {
		    if ( $safe_symbols ) {
			my $unimplemented = sub {
			    die "\$$k is read-only within String::Interpolate";
			};
			tie my $s, 'String::Interpolate::Func', {
			    FETCH => sub { "$$v" },
			    STORE => $unimplemented,
			};
			\$s;
		    } else {
			$v;
		    }
		} else {
		    if ( $safe_symbols ) {
			\ "$v";
		    } else {
			\$v;
		    }
		}
	    };
	}
    }


    unless ( $code ) {
	unless ( defined $string ) {
    	    croak("No string to interpolate");
	}

	$string = "BEGIN{import strict qw(refs subs); $$self->{pragmas}}; sub{<<$dlm\n$string\n$dlm\n}";

	if ( $safe ) {
	    no strict 'refs';
	    for ( 'String::Interpolate::Func::AUTOLOAD',
		  'warnings::unimport',
		  'warnings::import',
		  'strict::unimport',
		  'strict::import' ) {
		*{"${pkg}::$_"} = \&$_;
	    }
	    # Remove taint and generate a poor man's Safe::Hole
	    no warnings 'redefine';
	    *{"${pkg}::String::Interpolate::code"} = $safe->reval( $string =~ /(.*)/s );
	    $code = 1; # Just a flag in this case
	    # prevent_blessed_error_hack;
	} else {
	    $pkg ||= $$self->{defpgk};
	    $code = reval "package $pkg; $string";
	}
	if ( $@ ) {
	    return if $$self->{trap};
	    croak( $@ );
	}
	
	$$self->{code} = $code;
    };

    # Restore taint by appending null cut from $string
    if ( $safe ) {
	local $taint_flag = substr( $string, 0, 0 ); 
	local $safe_hole = $$self->{safe_hole};
	$string = $safe->reval('&String::Interpolate::code');
	# prevent_blessed_error_hack;
	if ( $@ ) {
	    return if $$self->{trap};
	    croak( $@ );
	}
    } else {
	$string = $$self->{trap} ? eval { &$code } : &$code; 
    }
    chop $string;

    # If we copied the lexicals then we must clean house to
    # avoid keeping them spuriously alive.
    $self->free_tmppkg if $$self->{lexicals};

    $string;
}

=back

=head2 Functional interface

For those heathens who don't like the OO interface.

=over 4

=item safe_interpolate

Exportable function equivalent to String::Interpolate->safe->exec(LIST).

=cut

sub safe_interpolate {
    __PACKAGE__->safe->exec(@_);
}

=item interpolate

Exportable function equivalent to
String::Interpolate->lexicals->exec(LIST).

=cut

sub interpolate {
    __PACKAGE__->lexicals->exec(@_);
}

=back

=head2 Ancillary methods

The following methods provide alternative interfaces and some fine
tuning capabilities.

=over 4

=item trap

Tells the String::Interpolate object whether or not to trap
exceptions.

    $i->trap;    # Enable trapping 
    $i->trap(1); # Enable trapping 
    $i->trap(0); # Disable trapping 

Returns the object so that it can be tagged on to constructor calls.

    my $i = String::Interpolate->safe->trap(0);

If the trap(0) method has not been called then trapping is enabled when
using a Safe compartment.

=cut

sub trap {
    my $self = shift;
    $self = $self->new unless ref $self;
    my $trap = shift;
    $$self->{trap} = defined $trap ? $trap : 1;
    $self;
}

=item unsafe_underscore

Tells the String::Interpolate object whether or not to use "unsafe
underscore" mode.  In this mode no precautions are taken to prevent
malicious code attempting to reach outside it's Safe compartment
through the $_ and %_ variables.

    $i->unsafe_underscore;    # Enable unsafe underscore mode
    $i->unsafe_underscore(1); # Enable unsafe underscore mode
    $i->unsafe_underscore(0); # Disable unsafe underscore mode

Returns the object so that it can be tagged on to constructor calls.

=cut

sub unsafe_underscore {
    my $self = shift;
    $self = $self->new unless ref $self;
    my $unsafe_underscore = shift;
    $$self->{unsafe_underscore} = defined $unsafe_underscore ? $unsafe_underscore : 1;
    $self;
}

=item unsafe_symbols

Tells the String::Interpolate object whether or not to use "unsafe
symbol" mode.  In this mode variables are simply shared with the Safe
compartment rather than being safely hidden behind variables tied to
blessed closures.  The setting of this flag as no effect when not
using a Safe compartment.

    $i->unsafe_symbols;    # Enable unsafe symbol mode
    $i->unsafe_symbols(1); # Enable unsafe symbol mode
    $i->unsafe_symbols(0); # Disable unsafe symbol mode

Returns the object so that it can be tagged on to constructor calls.

=cut

sub unsafe_symbols {
    my $self = shift;
    $self = $self->new unless ref $self;
    my $unsafe_symbols = shift;
    $$self->{unsafe_symbols} = defined $unsafe_symbols ? $unsafe_symbols : 1;
    $self;
}

=over 4

=item lexicals

This feature is EXPERIMENTAL.  Do not use it in real code.

Tells the String::Interpolate object whether or not to use the
PadWalker module to import all lexical variables from the calling
context into the temporary package or Safe compartment.  By default
this does not happen as it is conceptually ugly and quite expensive.

    $i->lexicals;     # Enable lexicals
    $i->lexicals(1)   # Enable lexicals 
    $i->lexicals(0);  # Disable lexicals

Returns the object so that it can be tagged on to constructor calls.

    my $i = String::Interpolate->safe->lexicals;

Enabling lexicals with a Safe compartment like this will give the code
read-only access to all your lexical variables.

Note that the lexicals used are those in scope at the final call that
performs the interpolation, not those in scope when the
String::Interpolate object is constructed.  Also you can't have your
cake and eat it.  If you cannot use this feature at the same time as
an explicit package or Safe compartment.

=cut

sub lexicals {
    my $self = shift;
    $self = $self->new unless ref $self;
    my $lexicals = shift;
    if ( ( $$self->{lexicals} = defined $lexicals ? $lexicals : 1 ) ) {
	delete $$self->{pkg};
	delete $$self->{safe};
    }
    $self;
}

=item package

Instructs the String::Interpolate object to forget its current Safe
compartment or namespace and use the specified one henceforth.  The
package name can be specified as a string, a GLOB or a GLOB reference.
The trailing '::' may be ommited.  With an undefined argument this
method instructs the object to use a new automatically allocated
temporary namespace.

The package method Returns the object so that it can be tagged on to
constructor calls.  It can also be used as a constructor.

    my $i = String::Interpolate->package('Q');   # Use namespace Q::
    $i->package;                                 # Use temporary namespace
    $i->package(*R);                             # Use namespace R::
    $i->package(\*S::);                          # Use namespace S::

Note that the last two forms are not commonly used as GLOB or GLOB
reference arguments passed to the exec(), new() or methods are
automatically passed on the the package() method.

=cut

sub package {
    my $self = shift;
    $self = $self->new unless ref $self;
    my $pkg = shift;
    $pkg = *$pkg if ref $pkg eq 'GLOB'; 
    ($pkg) = $pkg =~ /^\*?(?:main::(?!$))*(.*?)(?:::)?$/ or die;
    $self->free_tmppkg;
    delete $$self->{safe};
    delete $$self->{implicit_safe};
    delete $$self->{lexicals};
    $$self->{pkg} = $$self->{explicit_pkg} = $pkg;
    $self;
}

=item safe_hole

Tells the String::Interpolate object whether or not to use a
Safe::Hole object to wrap callbacks to subroutines specified in the
symbol mapping hash.  Without a Safe::Hole eval(), symbolic references
and method calls in callbacks won't function normally.

    my $i = String::Interpolate->safe->safe_hole;
    # Without a Safe::Hole Wibble::wobble() would be inaccessible
    $i->{FOO} = sub () { Wibble->wobble };

This feature only makes sense when evaluating in a Safe compartment
and you can only use it if you have the Safe::Hole module installed.

    $i->safe_hole;         # Enable use of Safe::Hole 
    $i->safe_hole(1);      # Enable use of Safe::Hole 
    $i->safe_hole(0);      # Disable use of Safe::Hole
    $i->safe_hole($hole);  # Use the Safe::Hole object $hole

This method can also be called implicitly as follows.

    $i->(\'SAFE HOLE');    # Enable use of Safe::Hole 
    $i->(\'NO_SAFE_HOLE'); # Disable use of Safe::Hole
    $i->($hole);           # Use the Safe::Hole object $hole

The safe_hole() method returns the object so that it can be tagged on
to constructor calls.

=cut

sub safe_hole {
    my $self = shift;
    $self = $self->new unless ref $self;
    my $safe_hole = shift;
    unless ( UNIVERSAL::isa( $safe_hole, 'Safe::Hole' )) {
	if ( $safe_hole || !defined $safe_hole ) {
	    unless ( eval { require Safe::Hole; 1 } ) {
		require Carp;
	        Carp::croak('String::Interpolate::safe_hole() requires Safe::Hole module');
	    }
	    $safe_hole = Safe::Hole->new(($Safe::Hole::VERSION > 0.09) ? ({}) : ());
	} else {
	    undef $safe_hole;
	}
    }	
    $$self->{safe_hole} = $safe_hole;
    $self;
}

=item pragma

Specify various options including Perl code to be complied in a
BEGIN{} block prior to compiling the string to be interpolated.  When
working in a Safe compartment, what you can do here is, of course,
highly limited.  In practice this is only useful for calling the
import() an unimport() methods on the warnings and strict modules.

For the most commonly used values, to control the handling of
interpolating undefined values, the following shorthands can also be
used:

  NOWARN => 'unimport warnings qw(uninitialized)'
  WARN   => ''
  FATAL  => 'import warnings FATAL => qw(uninitialized); import strict qw(vars)'

The default state for a newly created String::Interpolate object is
NOWARN.  All other warnings are enabled as are 'refs' and 'subs'
strictures.

You can call pragma() implicitly by passing SCALAR references to
exec().  Furthermore pragma('TRAP') is a synonym for trap(1) and
pragma('NO TRAP') is a synonym for trap(0).  Similarly for lexicals(),
unsafe_symbols(), unsafe_underscore() and safe_hole().  This makes the
following statements equivalent:

    $i->(\'FATAL',\'NO TRAP',\'SAFE SYMBOLS');
    $i->pragma('FATAL','NO_TRAP','NO UNSAFE_SYMBOLS');
    $i->pragma('FATAL')->trap(0)->unsafe_symbols(0);

The pragma() method returns the object so that it can be tagged on to
constructor calls.

=cut

sub pragma {
    my $self = shift;
    $self = $self->new unless ref $self;
    for my $pragma ( @_ ) {
	my ( $no, $method, $un) =
	    $pragma =~ /^(NO[ _]?)?(LEXICALS|TRAP|SAFE[_ ]HOLE|(?:((?:UN)?)SAFE[_ ](?:SYMBOLS|UNDERSCORE)))$/;
	if ( $method ) {
	    # For methods that start 'un' but for which the 'un' has been ommited
	    # reinstate the un and invert the sense of the 'no' prefix.
	    if ( defined $un && !$un ) {
		$no = !$no;
		$method = "UN$method";
	    }
	    $method =~ tr/ A-Z/_a-z/;
	    $self->$method(!$no + 0);
	} else {
	    $$self->{pragma} = $preset_pragma{$pragma} || $pragma;
	}
    }
    $self;
}

sub DESTROY {
    shift->free_tmppkg;
}

sub free_tmppkg {
    my $self = shift;
    delete $$self->{code};
    delete $$self->{safe} if $$self->{implicit_safe};
    if ( $$self->{tmppkg} ) {
	require Symbol;
        Symbol::delete_package( delete $$self->{tmppkg} );
    }
}

=item positionals 

Returns, as an lvalue, the reference to the array that holds the
values to use for the positional variables $1 and so on.

    my @p = qw ( one two three ); 
    my $i = new String::Interpolate \@p; 
    $i->positionals->[1] = "TWO";      # $p[1] = "TWO";
    $i->positionals = [ qw ( X Y ) ];  # Forget @p, use anon array
    undef $i->positionals;             # $1 etc. inherted from caller 

=cut

sub positionals : lvalue {
    my $self = shift;
    $$self->{pos};
}

sub ashash {
    my $self = shift;
    $self->exec({}) unless $$self->{map};
    $$self->{map}[-1];
}
 
package String::Interpolate::AsArray;

sub TIEARRAY { my ($class, $thing ) = @_; bless \$thing, $class }

sub STORE { ${${$_[0]}}->{pos}[$_[1]-1]=$_[2] }

sub FETCH { 
    require Carp; 
    Carp::croak('String::Interpolate objects STORE-only in ARRAY context');
}

*FETCHSIZE = \&FETCH;

# A private and very secretive class to give secure access to an object

package String::Interpolate::Func;

sub wrap_hash {
    my $class = shift;
    my ($k,$v) = @_;
    my $unimplemented = sub {
	die "%$k is read-only within String::Interpolate";
    };
    tie my %h, $class, {
	FETCH => sub { "$v->{$_[0]}" },
	STORE => $unimplemented,
	DELETE => $unimplemented,
	FIRSTKEY => sub { keys %$v; each %$v },
	NEXTKEY => sub { each %$v },
    };
    \%h;
}

sub TIEARRAY { 
    my $actions = $_[1];
    bless sub {	
	return unless my $action = $actions->{+shift};
	# Launder the argument list in case $action is wrapped by Safe::Hole
	# If the interpolated string was tainted then so are any arguments
	# passed from it.
	@_ = map { "$taint_flag$_" } @_;
	goto &$action unless $safe_hole;
	$safe_hole->call($action,@_);
    }, $_[0];
}

*TIEHASH = \&TIEARRAY;
*TIESCALAR = \&TIEARRAY;

sub AUTOLOAD { 
    my $self = shift;
    unshift @_ => our($AUTOLOAD) =~ /(\w+)$/;
    goto &$self; 
}

1;
__END__
=back
