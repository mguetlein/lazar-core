#!/usr/bin/perl

# known end of message indicator the server uses
$flag = "~EOT~\n";

# create socket
use IO::Socket;
my $sock = new IO::Socket::INET (
    PeerAddr => 'localhost',
    PeerPort => '30000',
    Proto => 'tcp'
);
die "Could not create socket: $!" unless $sock;

# send message
#while (<>) {

#print "\nSending message to daemon...";
#$msg = "c1ncc(cc1)N(N=O)C";
#print "['"; print $msg; print "'] --> daemon \n\n";
#print $sock $_;
#}
# read reply
#print "Receiving reply from daemon...\n";

$response = <$sock>;
while ($response ne "$flag") {
$rep.=$response;
$response = <$sock>;
}


print "['"; print $rep; print "'] <-- daemon \n\n";

#close socket
#print "Closing socket.i\n";
close($sock)
