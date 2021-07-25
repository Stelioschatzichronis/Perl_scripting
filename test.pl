 use IO::Handle;
 
# create a pipe for the communication parent -> child
pipe(FROM_PARENT, TO_CHILD)     or die "pipe: $!";
# create a second pipe for the communication child -> parent
pipe(FROM_CHILD,  TO_PARENT)    or die "pipe: $!";
 
# Set the "autoflush" option in both pipes, to make the information available
# for reading immediately after it has been written (disable buffering)
TO_CHILD->autoflush(1);
TO_PARENT->autoflush(1);
 
if ($pid = fork) {
    # This code is executed in the parent process
    # close the sides of the pipes not used in the parent
    close FROM_PARENT; close TO_PARENT;
 
    # Write a message to the child
    print TO_CHILD "Hi, I am the parent with process ID $$\n";
 
    # Read a message from the child
    chomp($line = <FROM_CHILD>);
    print "The parent process with PID $$ received a message: '$line'\n";
    close FROM_CHILD; close TO_CHILD;
    # Wait until the child exits
    waitpid($pid,0);
} else {
    # This code is executed in the child process.
    die "Error in fork: $!" unless defined $pid;
    # close the sides of the pipes not used by the child
    close FROM_CHILD; close TO_CHILD;
 
    # Read a message from the parent
    chomp($line = <FROM_PARENT>);
    print "The child process with PID $$ received a message: '$line'\n";
 
    # Write a message to the parent
    print TO_PARENT "Hi, this is the child process with PID $$\n";
    close FROM_PARENT; close TO_PARENT;
    exit;
}
