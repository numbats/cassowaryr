This package has been previously submitted with a dependence on the alphahull package. It was rejected because it had an ACM license, reflecting the tripack dependency in alphahull. Since then the alphahull package has been updated with the dependency on tripack replaced with interp package, and thus with a GPL licence. We have updated our package to depend on the new alphahull package and also now use a GPL license. This rectifies the issue with our original submission.

There is another a package for computing scagnostics available on CRAN, called scagnostics. However, this is built on C code that cannot be debugged, and rJava, which makes the scagnostics not readily available on all platforms. For this reason we have written the code to compute the scagnostics from scratch, which also will allow others to develop them further.


