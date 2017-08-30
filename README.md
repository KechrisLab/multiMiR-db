# multiMiR-db
Database and server code to support the multiMiR R package.

This contains two seperate components which include code(html,cgi)
for the web-server and R code to load data from the source files
into an SQL Database(MySQL).

The data from the 14 supported databases is not packaged with the
multiMiR R package due to its size.  Instead the R package submits an SQL
query to the server which performs the database query and returns
results. The code to accomplish this on the server is in the cgi 
folder.  The code to build database tables from the source files 
is in the database folder. The R package does allow you to change 
the URL of the server for requests so it is relatively easy to 
use the R package on your own server.

# Why might you want to use this code?
- Unsuported Organisms - you could modify the code to support 
additional organisms
- Improved performance - if the server is not performing 
adequately you can create a mirror that should in most cases 
perform better.  If you experience problems with our server performance
please let us know.

# Database code

R code to process the initial files available from each datasource and load into the MySQL database is 
provided.  You can run these individually or all at once.

* additional details

# Server code

The html is not nessecary for the package to function, 
but a configured webserver to support the CGI and also to serve 
cutoff files is required to setup your own server. If you need 
access to unsupported organisms or just want .  

## Steps to setup a server
  1. Install a Web Server(Apache recommended).
  1. Setup a database server(MySQL recommeded).
  1. Run database code to setup/load database tables.
  1. Setup and test the cgi and cutoff files.
  1. Setup R package and change the URL to point to your server.
  
