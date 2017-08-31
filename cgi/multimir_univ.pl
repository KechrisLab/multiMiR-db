#!C:/strawberry/perl/bin/perl.exe
use DBI;
use CGI ":all";


# ****************************************************************************************
# * Variables to edit for your server

$CGIurl="http://your.server.url/cgi-bin/";
# Current Version should always be updated to allow defaulting to current version.
# Less of an issue with current R package as it checks for versions and connects to latest 
# current version
my $dbName= "default_db_name";
my $dbHostname= "localhost";  #EDIT if not running on same machine
my $dbUser="user";
my $dbPassword="password";

# * End of server dependent variables
# *****************************************************************************************




my $query = param('query');		# the complete MySQL query
my $org = param('org');			# organism
my $mirna = param('mirna');		# miRNA name/ID
my $target = param('target');		# target gene name/ID
my $table = param('table');		# MySQL table to query
my $disease_drug = param('disease_drug');	# disease/drug name
my $p_cutoff = param('p_cutoff');		# cutoff for predicted targets
my $p_cutoff_type = param('p_cutoff_type');# cutoff type for predicted targets

print header;

#If version parameter is specified
if(length(param('dbName'))>0){
    $dbName=param('dbName');
}

# TODO
# If $query is provided, use it to query the database directly.
# If $query is not provided, build the query using other parameters
if ($query eq ""){

}


my @Col_names = [];
my $ref;
my $Nrow=0;
if($query ne ""){
    my $dbh = DBI->connect("DBI:mysql:".$dbName.":".$dbHostname,$dbUser,$dbPassword,{PrintError=>1,RaiseError=>1});
    # query the database
    my $sth = $dbh->prepare($query) or die "Couldn't prepare query: ".$query."\n". $dbh->errstr;
    $sth->execute() or die "Couldn't execute query: ".$query."\n" . $sth->errstr;
    @Col_names = @{ $sth->{NAME} };
    $ref = $sth->fetchall_arrayref;
    # save the summary table
    $Nrow = @{$ref};
    $sth->finish();
    $dbh->disconnect();
}


print "<html>\n\n<body>\n\n";
print "<div align=\"center\">\n";
print "<h2>multiMiR Summary</h2>\n";
print "<table align=\"center\" border=\" \">\n";
print "<tr align=\"center\" valign=\"top\">";
print "<th>Rows in SQL Result</th><th>SQL Query</th></tr>\n";
print "<tr align=\"center\" valign=\"top\">";
print "<td>$Nrow</td><td>$query</td></tr>\n";
print "</table>\n</div>\n\n";

print "<br><br>\n";

# print the result table
print "<div align=\"center\">\n";
print "<h2>multiMiR Result</h2>\n";
if ($Nrow > 0){
    print "<table align=\"center\" border=\" \">\n";
    print "<tr align=\"center\" valign=\"top\">";
    print "<th>$_</th>" for @Col_names;
    print "</tr>\n";
    for (my $i=0;$i<$Nrow;$i++){
	print "<tr align=\"center\" valign=\"top\">";
	print "<td>$_</td>" for @{$ref->[$i]};
	print "</tr>\n";
    }
    print "</table>\n";
}else {
    print "No result.<br>\n";
}
print "</div>\n\n";

print "</body>\n\n</html>\n";



