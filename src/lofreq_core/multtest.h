#ifndef MULTTEST_H
#define MULTTEST_H

/* -*- c-file-style: "k&r"; indent-tabs-mode: nil; -*-
 *
 * FIXME Copyright update
 *
 *********************************************************************/

typedef enum
{
     MTC_NONE,
     MTC_BONF,
     MTC_HOLMBONF,
     MTC_FDR
} mtc_type_t;

#define STR(name) # name

static char *mtc_type_str[] = {
    STR(MTC_NONE),
    STR(MTC_BONF),
    STR(MTC_HOLMBONF),
    STR(MTC_FDR),
};

int
mtc_str_to_type(char *t);


#endif
