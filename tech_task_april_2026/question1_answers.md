# Question 1: BASH

## Question 1.1: How to properly track errors in:

### Bash, when you call other Unix commands.
You can check the exit status of the command that has been run by
looking at the contents of `$?` -- if this is `0`, the command was
successful. If not, there was an error.
For more clues about the error, you can also direct errors logged
to `STDERR` to a file to enable examination of error logs following
an unsuccessful command execution, using `2>`, and/or watch this
during the run using `tail`.  It is also possible to redirect `STDERR`
to `STDOUT` to gather all logging in one stream, using `2>&1`.

### Python script:
Python allows `try`, `except`, `else`, and `finally` statements to enable
graceful handling of errors. The code in the `try` block is the code
you want to execute, and you can add handling of potential error
classes in any number of following `except` blocks. A final
`except` block can be added, with no specified error class, to
specify how to handle cases with unpredicted errors. The `else`
block following the `except` block is executed if there are no
errors, and a `finally` block can be added which is executed
whether or not there is an error. For example:

```
a = 2
try :
    b = int( input( "Enter a number: " )
    a = a/b
# If someone enters 0
except ZeroDivisionError :
    print( "Can't divide by 0!" )
# If someone enters a non-numeric value, int() won't work
except ValueError :
    print( "Please enter a number!" )
# If something else happens that we didn't think of
except :
    print( "There was an error" )
# If there weren't any errors
else :
    print( f"Division successful! Answer is {a}" )
finally :
    print( "Finished processing." )
```

### Perl script
There are a few ways to manage error handling in Perl. The basic built-in way
is to wrap code in an `eval` block and then check the contents of `$@` for
errors. There are also various modules available such as `Error.pm` and
`Try::Tiny`. `Try::Tiny` allows similar try/catch behaviour as above in Python.
For example:

```
#!/usr/bin/env perl

use strict;
use warnings;
use 5.10.0;

use Try::Tiny;

my $a = 2;

my $answer = try {
    say "Enter a number:";
    my $b = <STDIN>;
    $a = $a/$b;
} catch {
    warn "There was an error: $_";    # Error message sent to $_
} finally {
    say "Processing finished."
};
```

### Database commandline tool
When running `psql` in Bash, as with any other command, the exit status of the
command can be checked using `$?`. Within the SQL code, `psql` allows variables
such as `ON_ERROR_STOP` to halt a script and return non-zero exit status if
there's an error, and `ON_ERROR_ROLLBACK` to roll back to the status before the
failing statement was run and exit with non-zero status.

## Question 1.2


## Question 1.3
Comments are important in Bash because the syntax can be difficult to interpret
if you are unfamiliar with it, unlike a more human-readable language like
Python.
Comments are entered by prepending with a `#` character. For example:

```
# This is a special case which tells the OS which interpreter to use to parse the code.
#!/bin/bash     

# Start Jupyter virtual environment
source ~/jupyter/bin/activate

# Start Jupyter server
jupyter notebook
```
