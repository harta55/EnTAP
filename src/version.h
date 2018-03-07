//
// Created by harta on 2/25/18.
//

#ifndef ENTAP_VERSION_H
#define ENTAP_VERSION_H

#define LICENSE_YEAR      2017

// Comment this out if it is debug code
//#define RELEASE_BUILD

#define MAJOR_VERSION     0
#define MINOR_VERSION     8
#define BUILD_VERSION     0

#ifndef RELEASE_BUILD
    #define DEBUG_BUILD
#endif

#define TO_STR2(x)             #x
#define TO_STR(x)              TO_STR2(x)

#ifdef DEBUG_BUILD
    #define ENTAP_VERSION_STR       (TO_STR(MAJOR_VERSION) "." TO_STR(MINOR_VERSION) "." \
                                    TO_STR(BUILD_VERSION) "-DEBUG")
#else
    #define ENTAP_VERSION_STR      (TO_STR(MAJOR_VERSION) "." TO_STR(MINOR_VERSION) "." \
                                    TO_STR(BUILD_VERSION))
#endif

#endif //ENTAP_VERSION_H
