AC_INIT([<package>],[<version>])

# Find the compiler and compiler flags used by R.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_LANG(C++)
AC_PROG_CPP

# Search for a system library, in this case the zlib compression library (which
# contains the function 'deflate'). This sets the variable 'ac_cv_search_deflate'
# to what must be passed to ${PKG_LIBS}.
AC_SEARCH_LIBS(deflate, z, [], [AC_ERROR([The zlib library is required.])])
AC_CHECK_HEADERS(zlib.h, [], [AC_ERROR([The zlib library headers are required.])])

PROJ_CONFIG="pkg-config proj"

if `$PROJ_CONFIG --exists` ; then
  AC_MSG_NOTICE([pkg-config proj exists, will use it])
  proj_config_ok=yes
else
  proj_config_ok=no
fi

AC_ARG_WITH([proj-include],
    AC_HELP_STRING([--with-proj-include=DIR],
                   [location of proj header files]),
    [proj_include_path=$withval])
if test [ -n "$proj_include_path" ] ; then
   AC_SUBST([PROJ_CPPFLAGS],["-I${proj_include_path}"])
else
  if test "${proj_config_ok}" = yes; then
    PROJ_INCLUDE_PATH=`${PROJ_CONFIG} --cflags`
    AC_SUBST([PROJ_CPPFLAGS],["${PROJ_INCLUDE_PATH}"])
  fi
fi

# honor PKG_xx overrides
# for CPPFLAGS we will superfluously double R's flags
# since we'll set PKG_CPPFLAGS with this, but that shouldn't hurt

AC_ARG_WITH([proj-api],
    AC_HELP_STRING([--with-proj-api=yes/no],
               [use the deprecated proj_api.h even when PROJ 6 is available; default no]),
               [proj_api=$withval])
PROJ6="no"
PROJH="no"
if test "${proj_config_ok}" = yes; then
  PROJ_VERSION=`${PROJ_CONFIG} --modversion`
  PROJV1=`echo "${PROJ_VERSION}" | cut -c 1`
  if test "${PROJV1}" -ge 6; then
    PROJ6="yes"
    PROJ_CPPFLAGS="${PROJ_CPPFLAGS} -DHAVE_PROJ_H"
    if test "${proj_api}" = yes; then
        AC_MSG_NOTICE([using proj_api.h even with PROJ 5/6])
        PROJ_CPPFLAGS="${PROJ_CPPFLAGS} -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H"
    else
        AC_MSG_NOTICE([using proj.h.])
        PROJH="yes"
    fi
  fi
else
  if test "${PROJH}" = no ; then
    PROJH=yes
    AC_CHECK_HEADERS(proj.h,,PROJH=no)
    if test "${PROJH}" = yes; then
      PROJ6="yes"
      PROJ_CPPFLAGS="${PROJ_CPPFLAGS} -DHAVE_PROJ_H"
    fi
  fi
fi

CPPFLAGS="${INCPPFLAGS} ${PKG_CPPFLAGS} ${PROJ_CPPFLAGS}"

if test "${PROJH}" = no
then
  proj4ok=yes
  AC_CHECK_HEADERS(proj_api.h,,proj4ok=no)
  if test "${proj4ok}" = no; then
     AC_MSG_ERROR([proj_api.h not found in standard or given locations.])
  fi
fi

# dnl ditto for a library path
AC_ARG_WITH([proj-lib],
    AC_HELP_STRING([--with-proj-lib=LIB_PATH],
               [the location of proj libraries]),
               [proj_lib_path=$withval])
if test [ -n "$proj_lib_path" ] ; then
    AC_SUBST([PROJ_LIBS], ["-L${proj_lib_path} ${INPKG_LIBS} -lproj"])
else
  if test "${proj_config_ok}" = yes; then
    PROJ_LIB_PATH=`${PROJ_CONFIG} --libs`
    AC_SUBST([PROJ_LIBS], ["${PROJ_LIB_PATH} ${INPKG_LIBS}"])
	proj_version=`${PROJ_CONFIG} --modversion`
	AC_MSG_NOTICE([PROJ: ${proj_version}])
  else
    PROJ_LIBS="${PKG_LIBS} -lproj"
  fi
fi

LIBS="${PROJ_LIBS} ${INLIBS} ${PKG_LIBS}"

if test "${PROJH}" = no; then
  proj4ok=yes
  AC_CHECK_LIB(proj,pj_init_plus,,proj4ok=no)
  if test "${proj4ok}" = no; then
     AC_MSG_ERROR([libproj not found in standard or given locations.])
  fi
[cat > proj_conf_test.c <<_EOCONF
#include <stdio.h>
#include <stdlib.h>
#include <proj_api.h>
int main() {
    printf("%d\n", PJ_VERSION);
    exit(0);
}
_EOCONF]
else
[cat > proj_conf_test.cpp <<_EOCONF
#include <stdio.h>
#include <stdlib.h>
#include <proj.h>
int main() {
	proj_context_create();
    exit(0);
}
_EOCONF]
  #AC_CHECK_LIB(proj,proj_context_create,,proj6ok=no)
  AC_MSG_CHECKING(PROJ: checking whether PROJ and sqlite3 are available for linking:)
  ${CXX} ${CPPFLAGS} -o proj_conf_test proj_conf_test.cpp ${LIBS} -lsqlite3 2> errors.txt
  if test `echo $?` -ne 0 ; then
    proj6ok=no
    AC_MSG_RESULT(no)
  else
    proj6ok=yes
    AC_MSG_RESULT(yes)
  fi
  if test "${proj6ok}" = no; then
     AC_MSG_ERROR([libproj or sqlite3 not found in standard or given locations.])
  fi

[cat > proj_conf_test.c <<_EOCONF
#include <stdio.h>
#include <stdlib.h>
#include <proj.h>
int main() {
    printf("%d.%d.%d\n", PROJ_VERSION_MAJOR, PROJ_VERSION_MINOR, PROJ_VERSION_PATCH);
    exit(0);
}
_EOCONF]
fi

#AC_MSG_NOTICE([PKG_LIBS: ${PKG_LIBS}])

${CC} ${CFLAGS} ${CPPFLAGS} -o proj_conf_test proj_conf_test.c ${PROJ_LIBS}

proj_version=`./proj_conf_test`

AC_ARG_WITH([proj-share],
    AC_HELP_STRING([--with-proj-share=SHARE_PATH],
               [the location of proj metadata files]),
               [proj_share_path=$withval])
if test [ -n "$proj_share_path" ] ; then
    AC_MSG_NOTICE([PROJ_LIB: ${proj_share_path}])
fi

if test ${PROJ6} = "no"; then

[cat > proj_conf_test.c <<_EOCONF
#include <stdio.h>
#include <proj_api.h>
#if PJ_VERSION <= 480
FILE *pj_open_lib(projCtx, const char *, const char *);
#endif
int main() {
#if PJ_VERSION <= 480
    FILE *fp;
#else
    PAFile fp;
#endif
    projCtx ctx;
    ctx = pj_get_default_ctx();
    fp = pj_open_lib(ctx, "epsg", "rb");
    if (fp == NULL) exit(1);
#if PJ_VERSION <= 480
    fclose(fp);
#else
    pj_ctx_fclose(ctx, fp);
#endif
    exit(0);
}
_EOCONF]

${CC} ${CFLAGS} ${CPPFLAGS} -o proj_conf_test proj_conf_test.c ${PROJ_LIBS}
if test [ -n "$proj_share_path" ] ; then
  PROJ_LIB="${proj_share_path}" ./proj_conf_test
  proj_share=`echo $?`
else
  ./proj_conf_test
  proj_share=`echo $?`
fi

AC_MSG_CHECKING(PROJ: epsg found and readable)
if test ${proj_share} -eq 1 ; then
    AC_MSG_RESULT(no)
    STOP="stop"
else
    AC_MSG_RESULT(yes)
fi

rm -f proj_conf_test.c proj_conf_test

if test "$STOP" = "stop" ; then
    echo "Error: proj/epsg not found"
    echo "Either install missing proj support files, for example"
    echo "the proj-nad and proj-epsg RPMs on systems using RPMs,"
    echo "or if installed but not autodetected, set PROJ_LIB to the"
    echo "correct path, and if need be use the --with-proj-share="
    echo "configure argument."
    exit 1
fi

else # proj >= 6

if test "${PROJH}" = no; then

[cat > proj_conf_test.c <<_EOCONF
#include <stdio.h>
#include <proj_api.h>
int main() {
    PAFile fp;
    projCtx ctx;
    ctx = pj_get_default_ctx();
    fp = pj_open_lib(ctx, "proj.db", "rb");
    if (fp == NULL) exit(1);
    pj_ctx_fclose(ctx, fp);
    exit(0);
}
_EOCONF]

${CC} ${CFLAGS} ${CPPFLAGS} -o proj_conf_test proj_conf_test.c ${PROJ_LIBS}
if test [ -n "$proj_share_path" ] ; then
  PROJ_LIB="${proj_share_path}" ./proj_conf_test
  proj_share=`echo $?`
else
  ./proj_conf_test
  proj_share=`echo $?`
fi

AC_MSG_CHECKING(PROJ: proj.db found and readable)
if test ${proj_share} -eq 1 ; then
    AC_MSG_RESULT(no)
    STOP="stop"
else
    AC_MSG_RESULT(yes)
fi

rm -f proj_conf_test.c proj_conf_test

if test "$STOP" = "stop" ; then
    echo "Error: proj/proj.db not found"
    echo "Either install missing proj support files, set PROJ_LIB to the"
    echo "correct path, and if need be use the --with-proj-share="
    echo "configure argument."
    exit 1
fi

[cat > proj_conf_test.c <<_EOCONF
#include <stdio.h>
#include <proj_api.h>
#if PJ_VERSION <= 480
FILE *pj_open_lib(projCtx, const char *, const char *);
#endif
int main() {
#if PJ_VERSION <= 480
    FILE *fp;
#else
    PAFile fp;
#endif
    projCtx ctx;
    ctx = pj_get_default_ctx();
    fp = pj_open_lib(ctx, "conus", "rb");
    if (fp == NULL) exit(1);
#if PJ_VERSION <= 480
    fclose(fp);
#else
    pj_ctx_fclose(ctx, fp);
#endif
    exit(0);
}
_EOCONF]

${CC} ${CFLAGS} ${CPPFLAGS} -o proj_conf_test proj_conf_test.c ${PROJ_LIBS}
if test [ -n "$proj_share_path" ] ; then
  PROJ_LIB="${proj_share_path}" ./proj_conf_test
  proj_share=`echo $?`
else
  ./proj_conf_test
  proj_share=`echo $?`
fi

AC_MSG_CHECKING(PROJ: conus found and readable)
if test ${proj_share} -eq 1 ; then
    WARN="warn"
    AC_MSG_RESULT(no)
else
    AC_MSG_RESULT(yes)
fi

rm -f proj_conf_test.c proj_conf_test

if test "$WARN" = "warn" ; then
    echo "Note: proj/conus not found"
    echo "No support available in PROJ4 for NAD grid datum transformations"
    echo "If required, consider re-installing from source with the contents"
    echo "of proj-datumgrid-1.<latest>.zip from http://download.osgeo.org/proj/ in nad/."
fi

fi # PROJH = no

fi # proj >= 6


# Write the flags into the src/Makevars file.
AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS} ${ac_cv_search_deflate}"])
AC_SUBST([PKG_LIBS], ["${PROJ_LIBS} ${PKG_LIBS}"])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${PKG_CPPFLAGS}
    libs:     ${PKG_LIBS}

  --------------------------------------------------
"
