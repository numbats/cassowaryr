This is a new package submission

Note that there is a package for computing scagnostics available on CRAN, scagnostics. However, this is built on C code that cannot be debugged, and rJava, which makes the scagnostics not easily available on all platforms. For this reason we have written the code to compute them from scratch, which also will allow others to develop further.

We have addressed the problems with an initial submission. 

We have added a citation to a paper that describes some of the scagnostic methods. 
 
There was some confusion about licensing in the initial submission. We have changed the license to just read ACM. We notice that the alphahull package was re-published on CRAN Mar 26. Our package has the ACM license because it depends on alphahull. Ideally a new version of tripack will be available at some point in the future so that both packages can demove the ACM license dependency.  

