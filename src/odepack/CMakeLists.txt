enable_language(Fortran)

# There are some warnings since Odepack smartly uses some arrays as either
# integer or double arrays. We suppress those, since they are no error per se.
add_definitions(-std=legacy -w)

add_library(odepack opkdmain.f opkda1.f opkda2.f)
