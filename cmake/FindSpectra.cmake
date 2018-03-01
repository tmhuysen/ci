# FindSpectra.cmake. Find<package>.cmake-template from (https://cmake.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules).

# For the moment, only an installation of Spectra inside /usr/local is supported.

# Try to find Spectra
# Once done, this will define
#  spectra_FOUND            spectra is available on the system
#  spectra_INCLUDE_DIRS     the spectra include directories

# Since it's a header-only library, the include directories are enough.


# For the moment, only an installation in /usr/local is supported.
find_path(SPECTRA_PREFIX README.md HINTS /usr/local/spectra)

if("${SPECTRA_PREFIX}" STREQUAL "SPECTRA_PREFIX-NOTFOUND")
    message(WARNING "Spectra was not found in location /usr/local/.")
else()
    # If found, we let the user know that Spectra was found.
    message(STATUS "Spectra was found in /usr/local/spectra")

    # Set FOUND.
    set(spectra_FOUND TRUE)

    # Set the INCLUDE_DIRS.
    set(spectra_INCLUDE_DIRS "/usr/local/spectra/include")
endif()
