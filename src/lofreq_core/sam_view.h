/*
 * Modelled after original main_samview() variables and its getopts parsing
 * and added globals. then removed those not set in getops parsing
 */
typedef struct {
        int is_header;
        int is_header_only;
        int is_bamin;
        int compress_level;
        int is_bamout;
        int is_count;
        int of_type;
        int is_long_help;
        char *fn_out;
        char *fn_list;
        char *fn_ref;
        char *fn_rg;

        /* globals */
        int g_min_mapQ;
        int g_flag_on;
        int g_flag_off;
        float g_subsam;
        char *g_library;
        char *g_rg;
        /*void *g_bed;*/ char *bed_file;
} sam_view_opts_t;

int
new_sam_view_opts(sam_view_opts_t **svo);
void
free_sam_view_opts(sam_view_opts_t **svo);

#ifdef SAM_VIEW_ORIG
int main_samview(int argc, char *argv[]);
#else
int main_samview(char *bam_file, sam_view_opts_t *svo);
#endif
