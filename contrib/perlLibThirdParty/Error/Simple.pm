# Error.pm
#
# Copyright (c) 2006 Shlomi Fish <shlomif@iglu.org.il>. All rights reserved.
# This program is free software; you can redistribute it and/or
# modify it under the terms of the MIT/X11 license.
use Error;

1;
__END__
=head1 NAME

Error::Simple - the simple error sub-class of Error

=head1 SYNOPSIS

    use base 'Error::Simple';

=head1 DESCRIPTION

The only purpose of this module is to allow one to say:

    use base 'Error::Simple';

and the only thing it does is "use" Error.pm. Refer to the documentation
of L<Error> for more information about Error::Simple.

=head1 KNOWN BUGS

None.

=head1 AUTHORS

Shlomi Fish ( C<< shlomif@iglu.org.il >> )

=head1 SEE ALSO

L<Error>

