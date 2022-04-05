## Test environments
* local ubuntu 20.04 install, R 4.1.3
* rhubs Debian Linux, R-devel, clang, ISO-8859-15 locale (debian-clang-devel)
* rhubs Windows Server 2022, R-devel, 64 bit (windows-x86_64-devel)

## R CMD check results
There were no ERRORs, WARNINGs.

## Notes
I got an e-mail from Brian Ripley about an error in the current version on CRAN: https://www.stats.ox.ac.uk/pub/bdr/LENGTH1_self/accSDA.out

There are also a couple of notes in imports and lasydata.

All of these problems have been fixed in the current version.

Also removed links for badges from readme that are causing issues.
