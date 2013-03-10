library('mpeccable')

if(require("testthat", quietly = TRUE)) {
    test_package("mpeccable")
} else {
    print( "package 'testthat' not available, cannot run unit tests" )
}
