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
        
